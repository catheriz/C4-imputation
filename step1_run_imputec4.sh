#!/bin/bash

# Function to display usage
usage() {
    echo "Usage: $0 -i <input_vcf.gz> -o <output_prefix> [-b <build_version>] [-j <beagle_jar_path>]"
    echo "  -i  Input VCF file (required)"
    echo "  -o  Output prefix (required, will be used for output filenames)"
    echo "  -b  Genome build version (optional, default: 38, options: 37 or 38)"
    echo "  -j  Path to Beagle JAR file (optional, default: beagle.25Nov19.28d.jar)"
    echo "  -h  Show this help message"
    exit 1
}

# Default values
build=38
beagle_jar="beagle.25Nov19.28d.jar"

# Parse command-line arguments
while getopts "i:o:b:j:h" opt; do
    case ${opt} in
        i) vcf="$OPTARG" ;;
        o) out="$OPTARG" ;;
        b) build="$OPTARG" ;;
        j) beagle_jar="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Validate required arguments
if [[ -z "$vcf" || -z "$out" ]]; then
    echo "Error: Both input VCF (-i) and output prefix (-o) are required."
    usage
fi

# Validate build version
if [[ "$build" != "37" && "$build" != "38" ]]; then
    echo "Error: Invalid genome build. Use 37 or 38."
    exit 1
fi

# Check if the Beagle JAR file exists
if [[ ! -f "$beagle_jar" ]]; then
    echo "Error: Beagle JAR file not found at '$beagle_jar'. Please provide a valid path."
    exit 1
fi

# Define genomic region based on build version
declare -A reg=( ["37"]="chr6:24894177-33890574" ["38"]="chr6:24893949-33922797" )

echo "Processing input VCF: $vcf"
echo "Output prefix: $out"
echo "Using genome build: $build"
echo "Using Beagle JAR: $beagle_jar"

# Run bcftools and Beagle for imputation
bcftools view --no-version "$vcf" -r ${reg[$build]} | \
  java -Xmx8g -jar "$beagle_jar" \
  gt=/dev/stdin gp=true \
  ref=MHC_haplotypes_CEU_HapMap3_ref_panel.GRCh$build.vcf.gz \
  out="$out" \
  map=<(bcftools query -f "%CHROM\t%POS\n"MHC_haplotypes_CEU_HapMap3_ref_panel.GRCh$build.vcf.gz | \
  awk '{print $1"\t.\t"$2/1e7"\t"$2}')

echo "Processing complete. Output saved with prefix: $out"

