#!/bin/bash

# DESCRIPTION
#    Script for binning long reads based on UMIs. Part of 
#    longread_umi.
#    
# IMPLEMENTATION
#    author   Søren Karst (sorenkarst@gmail.com)
#             Ryan Ziels (ziels@mail.ubc.ca)
#    license  GNU General Public License
#
# TO DO
#    Add terminal messages.
#    Optimize trimming and filtering for speed.
#    Add bin size limit
#    Add mapping against adaptors to remove UMI artifacts

USAGE="
-- longread_umi read_binning_gawk: Longread UMI detection and read binning.
   Tool requires UMIs in both ends of the read flanked by defined
   adaptor regions. Perform just the gawk binning step.

usage: $(basename "$0" .sh) [-h] (-o dir)
( -p )
(-u value -U value -O value -S value) 

where:
    -h  Show this help text.
    -o  Output directory originally used in part 1 of umi binning.
    -p  Flag to disable Nanopore trimming and filtering.
        Use with PacBio reads.
    -u  Discard bins with a mean UMI match error above u.
    -U  Discard bins with a UMI match error standard
        deviation above U.
    -O  Normalize read orientation fraction to 'O' if < 'O' reads are
        either +/- strand orientation.
    -N  Max number of reads with +/- orientation. [Default = 10000]
    -S  UMI bin size/UMI cluster size cutoff. [Default = 10]
"

### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hzo:pu:U:O:N:S:' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    o) OUT_DIR=$OPTARG;;
    p) TRIM_FLAG=YES;;
    u) UMI_MATCH_ERROR=$OPTARG;;
    U) UMI_MATCH_ERROR_SD=$OPTARG;;
    O) RO_FRAC=$OPTARG;;
    N) MAX_BIN_SIZE=$OPTARG;;
    S) BIN_CLUSTER_RATIO=$OPTARG;;
    :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
    \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
  esac
done

# Check missing arguments
MISSING="is missing but required. Exiting."
if [ -z ${OUT_DIR+x} ]; then echo "-o $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${UMI_MATCH_ERROR+x} ]; then echo "-u $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${UMI_MATCH_ERROR_SD+x} ]; then echo "-U $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${RO_FRAC+x} ]; then echo "-O $MISSING"; echo ""; echo "$USAGE"; exit 1; fi;
if [ -z ${MAX_BIN_SIZE+x} ]; then echo "-N is missing. Defaulting to 10000 +/- reads ."; MAX_BIN_SIZE=10000; fi;
if [ -z ${BIN_CLUSTER_RATIO+x} ]; then echo "-S is missing. Defaulting to 10 ."; BIN_CLUSTER_RATIO=10; fi;

### Source commands and subscripts -------------------------------------
. $LONGREAD_UMI_PATH/scripts/dependencies.sh # Path to dependencies script

BINNING_DIR=$OUT_DIR/read_binning

# UMI binning and filtering

$GAWK \
  -v BD="$BINNING_DIR" \
  -v UME_MATCH_ERROR="$UMI_MATCH_ERROR" \
  -v UME_MATCH_ERROR_SD="$UMI_MATCH_ERROR_SD" \
  -v RO_FRAC="$RO_FRAC" \
  -v MAX_BIN_SIZE="$MAX_BIN_SIZE"  \
  -v BIN_CLUSTER_RATIO="$BIN_CLUSTER_RATIO" \
  '
  NR==1 {
    print "[" strftime("%T") "] ### Read-UMI match filtering ###" > "/dev/stderr";
    print "[" strftime("%T") "] Reading UMI1 match file..." > "/dev/stderr";
  }
  # Read UMI match file
  NR==FNR{
    # Extract data from optional fields
    for (i=12; i <= NF; i++){
      # Find NM field and remove prefix (primary hit err)
      if($i ~ /^NM:i:/){sub("NM:i:", "", $i); perr = $i};
      # Find secondary hit field, remove prefix and split hits
      if($i ~ /^XA:Z:/){sub("XA:Z:", "", $i); split($i, shits, ";")};
    }
    # Add primary mapping to hit list
    err1[$1][$3]=perr;
    # Add secondary mapping to hit list
    #Iterate over each hit
    for (i in shits){
      # Split hit in subdata (read, pos, cigar, err)  
      split(shits[i], tmp, ",");
      # Add hit if non-empty, not seen before and not target reverse strand
      if (tmp[1] != "" && !(tmp[1] in err1[$1]) && tmp[2] ~ "+"){
        err1[$1][tmp[1]] = tmp[4];
      }
    }
    next;
  }
  FNR==1 {
   print "[" strftime("%T") "] Reading UMI2 match file..." > "/dev/stderr";
  }
  # Read UMI match file
  {
    # Extract data from optional fields
    for (i=12; i <= NF; i++){
      # Find NM field and remove prefix (primary hit err)
      if($i ~ /^NM:i:/){sub("NM:i:", "", $i); perr = $i};
      # Find secondary hit field and remove prefix
      if($i ~ /^XA:Z:/){sub("XA:Z:", "", $i); split($i, shits, ";")};
    }
    # Add primary mapping to hit list
    err2[$1][$3]=perr;
    # Add secondary mapping to hit list
    # Split list of hits 
    #Iterate over each hit
    for (i in shits){
      # Split hit in subdata (read, pos, cigar, err)
      split(shits[i], tmp, ",");
      # Add hit if non-empty, not seen before and not target reverse strand
      if (tmp[1] != "" && !(tmp[1] in err2[$1]) && tmp[2] ~ "+"){
        err2[$1][tmp[1]] = tmp[4];
      }
    }
  } END {
    print "[" strftime("%T") "] UMI match filtering..." > "/dev/stderr"; 
    # Filter reads based on UMI match error
    for (umi in err1){    
      for (read in err1[umi]){
        # Define vars
        e1 = err1[umi][read];
        e2 = err2[umi][read];
        # Filter reads not matching both UMIs
        if (e1 != "" && e2 != ""){
          # Filter based on mapping error 
          if (e1 + e2 <= 6 && e1 <= 3 && e2 <= 3){
            # Add read to bin list or replace bin assignment if error is lower
            if (!(read in match_err)){
              match_umi[read] = umi;
              match_err[read] = e1 + e2;
            } else if (match_err[read] > e1 + e2 ){
              match_umi[read] = umi;
              match_err[read] = e1 + e2;
            } 
          }
        }
      }
    }
    print "[" strftime("%T") "] Read orientation filtering..." > "/dev/stderr";
    # Count +/- strand reads
    for (s in match_umi){
      UM=match_umi[s]
      sub("_rc", "", UM)
      # Read orientation stats
      ROC=match(match_umi[s], /_rc/)
      if (ROC != 0){
        umi_ro_plus[UM]++
        roc[s]="+"
      } else {
        umi_ro_neg[UM]++
        roc[s]="-"
      }
      # Count reads per UMI bin
      umi_n_raw[UM]++;
    }
    
    # Calculate read orientation fraction
    for (u in umi_ro_plus){
      # Check read orientation fraction
      if (umi_ro_plus[u] > 1 && umi_ro_neg[u] > 1){
        if (umi_ro_plus[u]/(umi_ro_neg[u]+umi_ro_plus[u]) < RO_FRAC ){
          rof_check[u]="rof_subset"
          rof_sub_neg_n[u] = umi_ro_plus[u]*(1/RO_FRAC-1)
          rof_sub_pos_n[u] = rof_sub_neg_n[u]
        } else if (umi_ro_neg[u]/(umi_ro_neg[u]+umi_ro_plus[u]) < RO_FRAC ){
          rof_check[u]="rof_subset"
          rof_sub_neg_n[u]=umi_ro_neg[u]*(1/RO_FRAC-1)
          rof_sub_pos_n[u]=rof_sub_neg_n[u]
        } else {
          rof_check[u]="rof_ok"
          rof_sub_neg_n[u]=MAX_BIN_SIZE
          rof_sub_pos_n[u]=MAX_BIN_SIZE
        }
      } else {
        rof_check[u]="rof_fail"
      }
    }
    
    # Subset reads
    for (s in match_umi){
      UMI_NAME=match_umi[s]
      sub("_rc", "", UMI_NAME)
      if(roc[s] == "+"){
        if(rof_sub_pos_n[UMI_NAME]-- > 0){
          ror_filt[s]=UMI_NAME
        }
      } else if (roc[s] == "-"){
        if(rof_sub_neg_n[UMI_NAME]-- > 0){
          ror_filt[s]=UMI_NAME
        }
      }
    }

    print "[" strftime("%T") "] UMI match error filtering..." > "/dev/stderr";
    # Calculate UME stats
    for (s in ror_filt){
      UM=ror_filt[s]
      # Count matching reads
      umi_n[UM]++;
      # UMI match error stats
      umi_me_sum[UM] += match_err[s]
      umi_me_sq[UM] += (match_err[s])^2
    }

    # Check UMI match error
    for (u in umi_n){
      UME_MEAN[u] = umi_me_sum[u]/umi_n[u]
      UME_SD[u] = sqrt((umi_me_sq[u]-umi_me_sum[u]^2/umi_n[u])/umi_n[u])
      if (UME_MEAN[u] > UME_MATCH_ERROR || UME_SD[u] > UME_MATCH_ERROR_SD){
        ume_check[u] = "ume_fail"
      } else {
        ume_check[u] = "ume_ok"
      }
    }

    print "[" strftime("%T") "] UMI bin/cluster size ratio filtering..." > "/dev/stderr";
    for (u in umi_n){
      CLUSTER_SIZE=u
      sub(".*;", "", CLUSTER_SIZE)
      bcr[u]=umi_n_raw[u]/CLUSTER_SIZE
      if (bcr[u] > BIN_CLUSTER_RATIO){
        bcr_check[u] = "bcr_fail"
      } else {
        bcr_check[u] = "bcr_ok"
      }
    }

    # Print filtering stats
    print "umi_name", "read_n_raw", "read_n_filt", "read_n_plus", "read_n_neg", \
      "read_max_plus", "read_max_neg", "read_orientation_ratio", "ror_filter", \
      "umi_match_error_mean", "umi_match_error_sd", "ume_filter", "bin_cluster_ratio", \
      "bcr_filter" > BD"/umi_binning_stats.txt"
    for (u in umi_n){
      print u, umi_n_raw[u], umi_n[u], umi_ro_plus[u], umi_ro_neg[u], \
        rof_sub_pos_n[u] + umi_ro_plus[u], rof_sub_neg_n[u] + umi_ro_neg[u], rof_check[u], \
        UME_MEAN[u], UME_SD[u], ume_check[u], bcr[u], bcr_check[u]\
        > BD"/umi_binning_stats.txt"
    }
	
    print "[" strftime("%T") "] Print UMI matches..." > "/dev/stderr"; 
    for (s in ror_filt){
      UMI_NAME=ror_filt[s]
      if( \
          ume_check[UMI_NAME] == "ume_ok" && \
          rof_check[UMI_NAME] == "rof_ok" && \
          bcr_check[UMI_NAME] == "bcr_ok" \
      ){print UMI_NAME, s, match_err[s]}
    }
    # Print to terminal
    print "[" strftime("%T") "] Done." > "/dev/stderr"; 
  }
' $BINNING_DIR/umi1_map.sam $BINNING_DIR/umi2_map.sam > $BINNING_DIR/umi_bin_map.txt