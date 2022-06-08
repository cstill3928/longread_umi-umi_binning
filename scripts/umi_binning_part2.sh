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
-- longread_umi umi_binning_part2: Longread UMI detection and read binning.
   Tool requires UMIs in both ends of the read flanked by defined
   adaptor regions. This part is meant to be used after UMI binning by gawk
   or python based method.

usage: $(basename "$0" .sh) [-h] (-d file -o dir -m value -M value )
(-s value -e value -f string -F string -r string -R string -p )
(-u value -U value -O value -S value -t value) 

where:
    -h  Show this help text.
    -o  Output directory.
    -f  Forward adaptor sequence. 
    -F  Forward primer sequence.
    -r  Reverse adaptor sequence.
    -R  Reverse primer sequence.
    -p  Flag to disable Nanopore trimming and filtering.
        Use with PacBio reads.
    -t  Number of threads to use.
"

### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hzo:f:F:r:R:pt:' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    o) OUT_DIR=$OPTARG;;
    f) FW1=$OPTARG;;
    F) FW2=$OPTARG;;
    r) RV1=$OPTARG;;
    R) RV2=$OPTARG;;
    p) TRIM_FLAG=YES;;
    t) THREADS=$OPTARG;;
    :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
    \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
  esac
done

# Check missing arguments
MISSING="is missing but required. Exiting."
if [ -z ${OUT_DIR+x} ]; then echo "-o $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${FW1+x} ]; then echo "-f $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${FW2+x} ]; then echo "-F $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${RV1+x} ]; then echo "-r $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${RV2+x} ]; then echo "-R $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${THREADS+x} ]; then echo "-t is missing. Defaulting to 1 thread."; THREADS=1; fi;


### Source commands and subscripts -------------------------------------
. $LONGREAD_UMI_PATH/scripts/dependencies.sh # Path to dependencies script


# Extract binned reads

umi_binning() {
  # Input
  local UMIMAP=$1
  local OUT=$2

  # Binning
  $GAWK -v out="$OUT" '
    BEGIN {g=1; outsub="./"out"/"g; system("mkdir \047" outsub "\047");}
    NR==FNR {
      # Get read name
      sub(";.*", "", $1);
      # Associate read name and umi match
      bin[$2]=$1;
      # Assign umi to a folder group if it has none
      if (foldergrp[$1] == ""){
        j++;
        if (j <= 4000){
          foldergrp[$1]=g;
        } else {
          j = 0;
          g++;
          foldergrp[$1]=g;
          outsub="./"out"/"g;
          system("mkdir \047" outsub "\047");
        }
      }
      next;
    }
    FNR%4==1 {
      read=substr($1,2);
      bin_tmp=bin[read]
      if ( bin_tmp != "" ){
        binfile=out"/"foldergrp[bin_tmp]"/"bin_tmp"bins.fastq";
        print > binfile;
        getline; print > binfile;
        getline; print > binfile;
        getline; print > binfile;
      }
    }
  ' $UMIMAP -
}

export -f umi_binning

#Make important directories
TRIM_DIR=$OUT_DIR/trim
BINNING_DIR=$OUT_DIR/read_binning

cat $TRIM_DIR/reads_tf.fq |\
  $GNUPARALLEL \
    --env umi_binning \
    -L4 \
	-j $THREADS \
	--block 300M \
	--pipe \
  "mkdir $BINNING_DIR/bins/job{#};\
  cat | umi_binning $BINNING_DIR/umi_bin_map.txt\
  $BINNING_DIR/bins/job{#}"

aggregate_bins() {
  # Input
  local IN=$1
  local OUTDIR=$2
  local OUTNAME=$3
  local JOB=$4

  # Determine output folder
  local BIN=$(( ($JOB - 1)/4000 ))
  mkdir -p $OUTDIR/$BIN

  # Aggregate data
  cat $IN > $OUTDIR/$BIN/$OUTNAME  
}

export -f aggregate_bins

find $BINNING_DIR/bins/*/*/ -name "*bins.fastq" -printf "%f\n" |\
  sort |\
  uniq |\
  $GNUPARALLEL \
    --env aggregate_bins \
    -j $THREADS \
	"aggregate_bins '$BINNING_DIR/bins/*/*/'{/} \
    $BINNING_DIR/bins {/} {#}"

rm -r $BINNING_DIR/bins/job*

## Testing
exit 0
# Filtering based on sub UMIs
cat \
  $UMI_DIR/reads_tf_start.fq \
  <($SEQTK seq -r $UMI_DIR/reads_tf_end.fq) |\
$CUTADAPT \
  -j $THREADS \
  -e 0.2 \
  -O 11 \
  -m 18 \
  -M 18 \
  --discard-untrimmed \
  -g $FW1...$FW2 \
  -g $RV1...$RV2 \
  - 2> $UMI_DIR/sumi_trim.log |\
  $GAWK '
    NR%4==1{print ">" substr($1, 2)}
    NR%4==2{print $0}
  ' > $UMI_DIR/sumi.fa

$USEARCH \
  -fastx_uniques $UMI_DIR/sumi.fa \
  -fastaout $UMI_DIR/sumi_u.fa \
  -relabel sumi \
  -strand plus \
  -sizeout \
  -minuniquesize 2

PATTERN="[ATCG]{3}[CT][AG][ATCG]{3}[CT][AG][ATCG]{3}[CT][AG][ATCG]{3}"
grep -B1 -E "$PATTERN" $UMI_DIR/sumi_u.fa |\
  sed '/^--$/d' > $UMI_DIR/sumi_f.fa

$SEQTK seq -r $UMI_DIR/sumi_f.fa > $UMI_DIR/sumi_frc.fa

$GAWK \
  -v UF="$UMI_DIR/sumi_f.fa" \
  -v UR="$UMI_DIR/sumi_frc.fa" \
  -v UFR="umi12u.fa" \
  '
  (FILENAME == UF && FNR%2==1){
    SIZE=$1
    sub(".*=|;", "", SIZE)
    getline
    UMI_FWD[$0]=SIZE
  }
  (FILENAME == UR && FNR%2==1){
    SIZE=$1
    gsub(".*=|;", "", SIZE)
    getline
    UMI_RV[$0]=SIZE
  }
  (FILENAME == UFR && FNR%2==1){
    NAME=$1
    getline
    # Extract UMI1 and UMI2 in both orientations
    s1 = substr($1, 1, 18);
    s2 = substr($1, 19, 36);
    if (s1 in UMI_FWD && s2 in UMI_RV){
      print NAME UMI_FWD[s1] ";" UMI_RV[s2] ";" s1 ";" s2 ";\n" $1 
    }
  }
 ' \
 $UMI_DIR/sumi_f.fa \
 $UMI_DIR/sumi_frc.fa \
 $UMI_DIR/umi12u.fa \
 > $UMI_DIR/umi12uf.fa