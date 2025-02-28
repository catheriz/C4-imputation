#!/bin/bash

# Load bcftools module
module load bcftools

# Function to display usage
usage() {
    echo "Usage: $0 -i <input_vcf.gz> [-o <output_vcf.gz>] [-d <dosage_file.txt>] [-h]"
    echo "  -i  Input VCF file (required)"
    echo "  -o  Output VCF file (optional, default: <input>_DS.vcf.gz)"
    echo "  -d  Dosage output file (optional, default: dosages.txt)"
    echo "  -h  Display this help message"
    exit 1
}

# Parse command-line arguments
while getopts "i:o:d:h" opt; do
    case ${opt} in
        i) INPUT_VCF="$OPTARG" ;;
        o) OUTPUT_VCF="$OPTARG" ;;
        d) DOSAGE_FILE="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Check if input VCF file is provided
if [ -z "$INPUT_VCF" ]; then
    echo "Error: Input VCF file is required."
    usage
fi

# Set default output filenames if not specified
if [ -z "$OUTPUT_VCF" ]; then
    OUTPUT_VCF="${INPUT_VCF%.vcf.gz}_DS.vcf.gz"
fi
if [ -z "$DOSAGE_FILE" ]; then
    DOSAGE_FILE="dosages.txt"
fi

# Annotate VCF file
echo "Processing: $INPUT_VCF"
echo "Output VCF file: $OUTPUT_VCF"
bcftools annotate -x ^INFO/DR2,^INFO/AF,^FORMAT/DS "$INPUT_VCF" -Oz -o "$OUTPUT_VCF"

# Extract ALT, DR2, AF, and dosage values
(echo -e "ALT\tDR2\tAF\t$(bcftools query -l "$OUTPUT_VCF" | tr '\n' '\t' | sed 's/\t$//')"; \
bcftools view -i 'ID=="C4"' "$OUTPUT_VCF" | bcftools query -f '%ALT\t%INFO/DR2\t%INFO/AF[\t%DS]\n') > "$DOSAGE_FILE"

echo "Processing complete. Output files:"
echo "  - VCF: $OUTPUT_VCF"
echo "  - Dosage file: $DOSAGE_FILE"
