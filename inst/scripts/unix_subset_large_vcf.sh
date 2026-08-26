#!/bin/bash

progname=$(basename "$0")

usage() {
cat << HEREDOC

 About:
   This script subsets a large VCF file. The largest chromosome will have 101
   sites in the subset VCF file. The number of sites for remaining chromosomes
   is proportional to their length.

 Usage:
   $progname -i INPUT_VCF_FILE -o OUTPUT_DIRECTORY

 Options:
   -i, --input-vcf     Full path to the initial VCF file
   -o, --output-dir    Full path where output files will be stored
   -h, --help          Display help message
   -v, --verbose       Verbose output

HEREDOC
}

inform() {
    TIMESTAMP=$(date +"%a %b %d %X %Y")
    echo "[$TIMESTAMP] $*" >&2
}

problem() {
    TIMESTAMP=$(date +"%a %b %d %T %Y")
    echo -e "\n[$TIMESTAMP] *ERROR*: $*" >&2
}

validate() {
    if [[ -z "$i" ]]; then
        problem "Input VCF (-i) not provided."
        exit 1
    fi
    if [[ ! -f "$i" ]]; then
        problem "Input VCF file not found: $i"
        exit 1
    fi
    inform "Input VCF: $i"

    if [[ -z "$o" ]]; then
        problem "Output directory (-o) not provided."
        exit 1
    fi
    if [[ ! -d "$o" ]]; then
        problem "Output directory not found: $o"
        exit 1
    fi
    inform "Output directory: $o"
}

# Initialize variables
verbose=0

# Parse options safely
OPTS=$(getopt -o "hi:o:v" --long "help,input-vcf:,output-dir:,verbose" -n "$progname" -- "$@")
if [[ $? != 0 ]]; then
    usage
    exit 1
fi

eval set -- "$OPTS"

while true; do
    case "$1" in
        -h|--help)
            usage
            exit 0
            ;;
        -i|--input-vcf)
            i="$2"
            shift 2
            ;;
        -o|--output-dir)
            o="$2"
            shift 2
            ;;
        -v|--verbose)
            verbose=$((verbose + 1))
            shift
            ;;
        --)
            shift
            break
            ;;
        *)
            break
            ;;
    esac
done

# Validate inputs
validate


# set a variable with the chromosome size
chrom_length=(640851 947102 1067971 1200490 1343557 1418242 1445207 1472805 1541735 1687656 2038340 2271494 2925236 3291936)

# Trouver longueur max
max_length=0
for len in "${chrom_length[@]}"; do
  (( len > max_length )) && max_length=$len
done

# Fonction ceiling_divide qui arrondit à l'entier supérieur quotient*100/max_length
ceiling_divide() {
  dividend=$1
  divisor=$2
  quotient=$(( dividend / divisor ))
  remainder=$(( dividend % divisor ))
  if (( remainder > 0 )); then
    (( quotient += 1 ))
  fi
  echo $quotient
}

# Calculate the subset sizes
subsets=()
for len in "${chrom_length[@]}"; do
  scaled=$(( len * 100 ))
  subset=$(ceiling_divide $scaled $max_length)
  subsets+=($subset)
done

# set a variable with the chromosome names
chroms=("Pf3D7_01_v3" "Pf3D7_02_v3" "Pf3D7_03_v3" "Pf3D7_04_v3" \
        "Pf3D7_05_v3" "Pf3D7_06_v3" "Pf3D7_07_v3" "Pf3D7_08_v3" \
        "Pf3D7_09_v3" "Pf3D7_10_v3" "Pf3D7_11_v3" "Pf3D7_12_v3" \
        "Pf3D7_13_v3" "Pf3D7_14_v3")

# set the path to your VCF file
# vcf_file="/Users/karimmane/Documents/Karim/Personnel/Transmission_paper/mpbr/inst/extdata/Input_Data.vcf.gz"
        
for (( i=0; i <= ${#subsets[@]}; i++ )); do
  echo ${subsets[i]} 
  j=$((i + 1))
  # Extract variants for this chromosome to temporary VCF
  bcftools view -r ${chroms[i]} $i -O z -o $o/chrom_$j.vcf.gz 
  tabix -p vcf $o/chrom_$j.vcf.gz
  
  # Extract variant positions or IDs
  chrom_pos=chrom_${j}_positions.txt
  bcftools query -r ${chroms[i]} -f '%POS\n' $o/chrom_$j.vcf.gz > $o/$chrom_pos

  # Shuffle positions and take first subset_size positions
  sampled_chrom_pos=chrom_${j}_positions_sampled.txt
  shuf $o/$chrom_pos | head -n ${subsets[i]} | sort -n > $o/$sampled_chrom_pos

  # Create a bed file to store the subset positions
  bed_file=chrom_${j}_regions.bed
  awk '{pos=$1; print "'${chroms[i]}'\t"(pos-1)"\t"pos; print "'${chroms[i]}'\t"pos"\t"(pos+1)}' $o/$sampled_chrom_pos > $o/$bed_file

  # Extract variants at sampled positions into final subset VCF for this chromosome
  bcftools view -R $o/$bed_file $o/chrom_$j.vcf.gz -O z -o $o/chrom_${j}_sampled.vcf.gz
done


# concatenate and index the subset VCF files
bcftools concat $o/*_sampled.vcf.gz -Oz -o $o/test_data.vcf.gz
tabix -p vcf $o/test_data.vcf.gz



