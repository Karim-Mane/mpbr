#!/bin/bash

progname=$(basename "$0")

usage() {
    cat << HEREDOC

 About:
   This script is designed for subsetting a large VCF file. The largest chromosome
   will have 101 sites in the subset VCF file. The number of sites for the remaining
   chromosomes is proportional to their length.

 Usage:
   $progname -v VCF_FILE -o OUTPUT_DIRECTORY

 Options:
   -v, --vcf           Full path to the initial VCF file
   -o, --output-dir    Full path to the folder where the output files will be stored
   -h, --help          Display this help message and exit
   -v, --verbose       Verbose

HEREDOC
}

inform() {
  TIMESTAMP=$(date +"%a %b %d %X %Y")
  echo "[$TIMESTAMP] $*" 1>&2
}

problem() {
  TIMESTAMP=$(date +"%a %b %d %T %Y")
  echo -e "\n[$TIMESTAMP] *ERROR*: $*" 1>&2
}

validate() {
    if [[ -z "$v" ]]; then
        problem "Input VCF file (-v/--vcf) is required."
        exit 1
    fi

    if [[ ! -f "$v" ]]; then
        problem "Input VCF file not found: $v"
        exit 1
    fi
    inform "The input VCF file is: $v"

    if [[ -z "$o" ]]; then
        problem "Output directory (-o/--output-dir) is required."
        exit 1
    fi

    if [[ ! -d "$o" ]]; then
        problem "Output directory not found: $o"
        exit 1
    fi
    inform "The output directory is: $o"
}

# ---------- main ----------

# No args?
if [[ $# -eq 0 ]]; then
    echo "No arguments were passed! See usage below."
    usage
    exit 1
fi

#verbose=0
v=""
o=""

# Manual argument parsing (portable)
while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            usage
            exit 0
            ;;
        -v|--vcf)
            if [[ -z "$2" ]]; then
                problem "Option $1 requires an argument."
                exit 1
            fi
            v="$2"
            shift 2
            ;;
        -o|--output-dir)
            if [[ -z "$2" ]]; then
                problem "Option $1 requires an argument."
                exit 1
            fi
            o="$2"
            shift 2
            ;;
        --) # explicit end of options
            shift
            break
            ;;
        -*)
            problem "Unknown option: $1"
            usage
            exit 1
            ;;
        *)
            # first non-option argument – stop parsing options
            break
            ;;
    esac
done

# Validate required inputs
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

# #for (( i=0; i < 1; i++ )); do        
for (( i=0; i < ${#subsets[@]}; i++ )); do 
  j=$((i + 1))
  # Extract variants for this chromosome to temporary VCF
  bcftools view -r ${chroms[i]} $v -O z -o $o/chrom_$j.vcf.gz 
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