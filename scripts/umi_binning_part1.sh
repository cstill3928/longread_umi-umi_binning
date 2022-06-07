#!/bin/bash

# DESCRIPTION
#    Script for binning long reads based on UMIs. Part of 
#    longread_umi.
#    
# IMPLEMENTATION
#    author   Søren Karst (sorenkarst@gmail.com)
#             Ryan Ziels (ziels@mail.ubc.ca)
#    modifier Chris Still (chris.still@epic-bio.com)
#    license  GNU General Public License
#
# TO DO
#    Add terminal messages.
#    Optimize trimming and filtering for speed.
#    Add bin size limit
#    Add mapping against adaptors to remove UMI artifacts

USAGE="
-- longread_umi umi_binning_part1: Longread UMI detection and read binning.
   Tool requires UMIs in both ends of the read flanked by defined
   adaptor regions. This part is meant to prepare the sam files used in umi binning 
   by the gawk or python based methods.

usage: $(basename "$0" .sh) [-h] (-d file -o dir -m value -M value )
(-s value -e value -f string -F string -r string -R string -p )
(-u value -U value -O value -S value -X value -t value) 

where:
    -h  Show this help text.
    -d  Reads in fastq format.
    -o  Output directory.
    -m  Minimum read length.
    -M  Maximum read length.
    -s  Check start of read up to s bp for UMIs.
    -e  Check end of read up to f bp for UMIs.
    -f  Forward adaptor sequence. 
    -F  Forward primer sequence.
    -r  Reverse adaptor sequence.
    -R  Reverse primer sequence.
    -p  Flag to disable Nanopore trimming and filtering.
        Use with PacBio reads.
    -t  Number of threads to use.
    -X  Number of parallel jobs for splitting up reads for porechop
"

### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hzd:o:m:M:s:e:f:F:r:R:pt:X:' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    d) READ_IN=$OPTARG;;
    o) OUT_DIR=$OPTARG;;
    m) MIN_LENGTH=$OPTARG;;
    M) MAX_LENGTH=$OPTARG;;
    s) START_READ_CHECK=$OPTARG;;
    e) END_READ_CHECK=$OPTARG;;
    f) FW1=$OPTARG;;
    F) FW2=$OPTARG;;
    r) RV1=$OPTARG;;
    R) RV2=$OPTARG;;
    p) TRIM_FLAG=YES;;
    X) PORECHOP_GNU_SPLIT=$OPTARG;;
    t) THREADS=$OPTARG;;
    :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
    \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
  esac
done

# Check missing arguments
MISSING="is missing but required. Exiting."
if [ -z ${READ_IN+x} ]; then echo "-d $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${OUT_DIR+x} ]; then echo "-o $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${MIN_LENGTH+x} ]; then echo "-m $MISSING"; echo ""; echo "$USAGE"; exit 1; fi; 
if [ -z ${MAX_LENGTH+x} ]; then echo "-M $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${START_READ_CHECK+x} ]; then echo "-s $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${END_READ_CHECK+x} ]; then echo "-e $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${FW1+x} ]; then echo "-f $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${FW2+x} ]; then echo "-F $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${RV1+x} ]; then echo "-r $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${RV2+x} ]; then echo "-R $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${PORECHOP_GNU_SPLIT+x} ]; then echo "-X is missing. Defaulting to 30 ."; PORECHOP_GNU_SPLIT=30; fi;
if [ -z ${THREADS+x} ]; then echo "-t is missing. Defaulting to 1 thread."; THREADS=1; fi;


### Source commands and subscripts -------------------------------------
. $LONGREAD_UMI_PATH/scripts/dependencies.sh # Path to dependencies script

### Primer formating
revcom() {
  echo $1 |\
  $GAWK '{print ">dummy\n" $0}' |\
  $SEQTK seq -r - |\
  $GAWK '!/^>/'  
}
FW1R=$(revcom "$FW1")
FW2R=$(revcom "$FW2")
RV1R=$(revcom "$RV1")
RV2R=$(revcom "$RV2")

### Read trimming and filtering -----------------------------------------------
mkdir $OUT_DIR
TRIM_DIR=$OUT_DIR/trim
mkdir $TRIM_DIR

# Trim data
if [ -z ${TRIM_FLAG+x} ]; then

  # Perform porechop and filtlong in parallel
  FT_THREADS=$(( $THREADS/$PORECHOP_GNU_SPLIT ))
  if (( FT_THREADS < 1 )); then
    FT_THREADS=1
  elif (( FT_THREADS > THREADS )); then
    FT_THREADS=1
  fi

  cat $READ_IN | $GNUPARALLEL --progress -j $PORECHOP_GNU_SPLIT -L 4 --round-robin --pipe \
    "cat > $TRIM_DIR/{#}.tmp;\
    $PORECHOP_UMI \
      -i $TRIM_DIR/{#}.tmp \
      -o $TRIM_DIR/{#}_trim.tmp \
      --threads $FT_THREADS \
      --min_split_read_size $MIN_LENGTH \
      --adapter_threshold  80 \
	  --min_trim_size 20 \
      --extra_end_trim 0 \
      --extra_middle_trim_good_side 0 \
      --extra_middle_trim_bad_side 0 \
      --middle_threshold 80 \
      --check_reads 5000; \
    $FILTLONG --min_length $MIN_LENGTH --min_mean_q 70 $TRIM_DIR/{#}_trim.tmp |\
      $CUTADAPT -j $FT_THREADS -m $MIN_LENGTH -M $MAX_LENGTH - \
        -o $TRIM_DIR/{#}_filt.tmp;"

  # Concatenate temp files
  cat $TRIM_DIR/*_filt.tmp > $TRIM_DIR/reads_tf.fq
  rm $TRIM_DIR/*.tmp
else
# Create symlink if already trimmed.
  ln -s $PWD/$READ_IN $PWD/$TRIM_DIR/reads_tf.fq  
fi

### Extract UMI references sequences ------------------------------------------- 
mkdir $OUT_DIR/umi_ref
UMI_DIR=$OUT_DIR/umi_ref

# Extract UMI terminal region
$GAWK -v UD="$UMI_DIR" 'NR%4==1{
       print $0 > UD"/reads_tf_start.fq";
       print $0 > UD"/reads_tf_end.fq";  
     }
     NR%4==2{
       print substr($0, 1, 200) > UD"/reads_tf_start.fq";
       print substr($0, length($0) - 199, 200)  > UD"/reads_tf_end.fq";  
     }
     NR%4==3{
       print $0 > UD"/reads_tf_start.fq";
       print $0 > UD"/reads_tf_end.fq";   
     }
     NR%4==0{
       print substr($0, 1, 200) > UD"/reads_tf_start.fq";
       print substr($0, length($0) - 199, 200)  > UD"/reads_tf_end.fq";  
     }
' $TRIM_DIR/reads_tf.fq


# Extract UMI pairs with correct lengths
$CUTADAPT -j $THREADS -e 0.2 -O 11 -m 18 -M 18 \
  --discard-untrimmed \
  -g $FW1...$FW2 -g $RV1...$RV2 \
  -G $RV2R...$RV1R -G $FW2R...$FW1R \
  -o $UMI_DIR/umi1.fq -p $UMI_DIR/umi2.fq \
  $UMI_DIR/reads_tf_start.fq $UMI_DIR/reads_tf_end.fq \
  > $UMI_DIR/perfect_trim.log

paste -d "" <( sed -n '1~4s/^@/>/p;2~4p' $UMI_DIR/umi1.fq ) \
            <( sed -n '1~4s/^@/>/p;2~4p' $UMI_DIR/umi2.fq ) |\
  cut -d " " -f1 > $UMI_DIR/umi12.fa

# Extract UMI pairs with correct patterns 

# Pattern: (NNNYRNNNYRNNNYRNNN NNNYRNNNYRNNNYRNNN)
PATTERN="[ATCG]{3}[CT][AG][ATCG]{3}[CT][AG][ATCG]{3}[CT][AG][ATCG]{6}\
[CT][AG][ATCG]{3}[CT][AG][ATCG]{3}[CT][AG][ATCG]{3}"
grep -B1 -E "$PATTERN" $UMI_DIR/umi12.fa |\
  sed '/^--$/d' > $UMI_DIR/umi12f.fa

# Cluster UMI pairs
$USEARCH \
  -fastx_uniques \
  $UMI_DIR/umi12f.fa \
  -fastaout $UMI_DIR/umi12u.fa \
  -sizeout \
  -minuniquesize 1 \
  -relabel umi \
  -strand both

$USEARCH \
  -cluster_fast $UMI_DIR/umi12u.fa \
  -id 0.90 \
  -centroids $UMI_DIR/umi12c.fa \
  -uc $UMI_DIR/umi12c.txt \
  -sizein \
  -sizeout \
  -strand both \
  -minsize 1

# Extract putative UMI pairs
$CUTADAPT -j $THREADS -e 0.2 -O 11 -m 18 -l 18 \
  --discard-untrimmed \
  -g $FW1 -g $RV1 \
  -G $RV2R -G $FW2R \
  -o $UMI_DIR/umi1p.fq -p $UMI_DIR/umi2p.fq \
  $UMI_DIR/reads_tf_start.fq $UMI_DIR/reads_tf_end.fq \
  > $UMI_DIR/putative_trim.log

paste -d "" <( sed -n '1~4s/^@/>/p;2~4p' $UMI_DIR/umi1p.fq ) \
            <( sed -n '1~4s/^@/>/p;2~4p' $UMI_DIR/umi2p.fq ) |\
  cut -d " " -f1 > $UMI_DIR/umi12p.fa

$BWA index $UMI_DIR/umi12c.fa

$BWA aln \
  $UMI_DIR/umi12c.fa \
  $UMI_DIR/umi12p.fa \
  -n 6 \
  -t $THREADS \
  -N > $UMI_DIR/umi12p_map.sai
$BWA samse \
  -n 10000000 \
  $UMI_DIR/umi12c.fa \
  $UMI_DIR/umi12p_map.sai \
  $UMI_DIR/umi12p.fa|\
  $SAMTOOLS view -F 4 - \
  > $UMI_DIR/umi12p_map.sam

$GAWK \
  -v UMS="$UMI_DIR/umi12p_map.sam" \
  -v UC="$UMI_DIR/umi12c.fa" \
  '
  (FILENAME == UMS){
      CLUSTER[$3]++
  }
  (FILENAME == UC && FNR%2==1){
    NAME=substr($1,2)
    if (NAME in CLUSTER){
      if (CLUSTER[NAME]+0 > 2){
        SIZE=CLUSTER[NAME]
        gsub(";.*", "", NAME)
        print ">" NAME ";size=" SIZE ";"
        getline; print
      }
    }
  }
  ' \
  $UMI_DIR/umi12p_map.sam \
  $UMI_DIR/umi12c.fa \
  > $UMI_DIR/umi12cf.fa 

# Remove potential chimeras
paste <(cat $UMI_DIR/umi12cf.fa | paste - - ) \
  <($GAWK '!/^>/{print}' $UMI_DIR/umi12cf.fa | rev | tr ATCG TAGC) |\
  $GAWK -v UD="$UMI_DIR" 'NR==FNR {
      #Format columns
      split($1, a, /[>;]/);
      sub("size=", "", a[3]);
      # Extract UMI1 and UMI2 in both orientations
      s1 = substr($2, 1, 18);
      s2 = substr($2, 19, 36);
      s1rc= substr($3, 1, 18);
      s2rc= substr($3, 19, 36);
      # Register UMI1 size if larger than current highest or if empty
      if ((g1n[s1]+0) <= (a[3]+0) || g1n[s1] == ""){
        g1n[s1] = a[3];
        g1[s1] = a[2];
      }
      # Register UMI2 size if larger than current highest or if empty
      if ((g2n[s2]+0) <= (a[3]+0) || g2n[s2] == ""){
        g2n[s2] = a[3];
        g2[s2] = a[2];
      }
      # Register UMI1rc size if larger than current highest or if empty
      if ((g1n[s1rc]+0) <= (a[3]+0) || g1n[s1rc] == ""){
        g1n[s1rc] = a[3];
        g1[s1rc] = a[2];
      }
      # Register UMI2rc size if larger than current highest or if empty
      if ((g2n[s2rc]+0) <= (a[3]+0) || g2n[s2rc] == ""){
        g2n[s2rc] = a[3];
        g2[s2rc] = a[2];
      }
      # Register UMI1 and UMI matches for current UMI
      u[a[2]] = a[3];
      s1a[a[2]] = s1;
      s2a[a[2]] = s2;
      s1arc[a[2]] = s1rc;
      s2arc[a[2]] = s2rc;
    } END {
      for (i in u){
        keep="no";
        if (g1[s1a[i]] == i && g2[s2a[i]] == i && g1[s1arc[i]] == i && g2[s2arc[i]] == i && s1a[i] != s1arc[i]){
          keep="yes";
          print ">"i";"u[i]"\n"s1a[i]s2a[i] > UD"/umi_ref.fa";
        } else if (s1a[i] == s1arc[i]){
          keep="tandem"
          print ">"i";"u[i]"\n"s1a[i]s2a[i] > UD"/umi_ref.fa";
        }
        print i, n[i], s1a[i], s2a[i], keep, g1[s1a[i]]"/"g2[s2a[i]]"/"g1[s1arc[i]]"/"g2[s2arc[i]], u[i]
      }  
    }' > $UMI_DIR/umi_ref.txt

### Bin reads based on UMIs ----------------------------------------------------
mkdir $OUT_DIR/read_binning
mkdir $OUT_DIR/read_binning/bins
BINNING_DIR=$OUT_DIR/read_binning

# Extract UMI region
$GAWK -v BD="$BINNING_DIR" -v TL="$START_READ_CHECK" '
  NR%4==1{
    print ">" substr($1,2) > BD"/reads_tf_umi1.fa";
  }
  NR%4==2{
    print substr($0, 1, TL) > BD"/reads_tf_umi1.fa";
  }
' $UMI_DIR/reads_tf_start.fq

$GAWK -v BD="$BINNING_DIR" -v TL="$END_READ_CHECK" '
  NR%4==1{
     print ">" substr($1,2) > BD"/reads_tf_umi2.fa";  
   }
   NR%4==2{
     print substr($0, length($0) - TL + 1, TL)  > BD"/reads_tf_umi2.fa";  
   }
' $UMI_DIR/reads_tf_end.fq


# Divide in barcode1 and barcode2 files
cat $UMI_DIR/umi_ref.fa <($SEQTK seq -r $UMI_DIR/umi_ref.fa |\
  $GAWK 'NR%2==1{print $0 "_rc"; getline; print};') |\
  $GAWK -v BD="$BINNING_DIR" 'NR%2==1{
       print $0 > BD"/umi_ref_b1.fa";
       print $0 > BD"/umi_ref_b2.fa";  
     }
     NR%2==0{
       print substr($0, 1, 18) > BD"/umi_ref_b1.fa";
       print substr($0, 19, 18)  > BD"/umi_ref_b2.fa";  
     }'

# Map UMIs to UMI references
## Important settings:
## -N : diasble iterative search. All possible hits are found.
## -F 20 : Removes unmapped and reverse read matches. Keeps UMIs
##         in correct orientations.

$BWA index $BINNING_DIR/reads_tf_umi1.fa
$BWA index $BINNING_DIR/reads_tf_umi2.fa

$BWA aln $BINNING_DIR/reads_tf_umi1.fa $BINNING_DIR/umi_ref_b1.fa \
  -n 3 -t $THREADS -N > $BINNING_DIR/umi1_map.sai
$BWA samse -n 10000000 $BINNING_DIR/reads_tf_umi1.fa $BINNING_DIR/umi1_map.sai\
  $BINNING_DIR/umi_ref_b1.fa | $SAMTOOLS view -F 20 - > $BINNING_DIR/umi1_map.sam

$BWA aln $BINNING_DIR/reads_tf_umi2.fa $BINNING_DIR/umi_ref_b2.fa \
  -n 3 -t $THREADS -N > $BINNING_DIR/umi2_map.sai
$BWA samse -n 10000000 $BINNING_DIR/reads_tf_umi2.fa $BINNING_DIR/umi2_map.sai\
  $BINNING_DIR/umi_ref_b2.fa | $SAMTOOLS view -F 20 - > $BINNING_DIR/umi2_map.sam