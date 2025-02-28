#!/usr/bin/env python
# coding: utf-8

import os
import argparse
import pandas as pd

# Argument parser to allow user input
parser = argparse.ArgumentParser(description="Process C4 haplotype dosages from an input file.")
parser.add_argument("-i", "--input", required=True, help="Input file name after step3")
parser.add_argument("-o", "--output_prefix", required=True, help="Output file prefix")
args = parser.parse_args()

# Check if the input file exists
if not os.path.exists(args.input):
    print(f"Error: The file '{args.input}' does not exist.")
    exit(1)
# Ensure output directory exists
output_dir = os.path.dirname(args.output_prefix)
if output_dir and not os.path.exists(output_dir):
    print(f"Error: The output directory '{output_dir}' does not exist.")
    exit(1)
# Generate output filenames
counts_file = f"{args.output_prefix}_counts.csv"
dosages_file = f"{args.output_prefix}_dosages.csv"

# Load the dosages data from user-specified file
dosages_info = pd.read_csv(args.input, sep='\t')

haplotypes = {
    '<AL-AL-1>': {'C4A': 2, 'C4B': 0, 'C4L': 2, 'C4S': 0},
    '<AL-AL-2>': {'C4A': 2, 'C4B': 0, 'C4L': 2, 'C4S': 0},
    '<AL-AL-3>': {'C4A': 2, 'C4B': 0, 'C4L': 2, 'C4S': 0},
    '<AL-BS-1>': {'C4A': 1, 'C4B':1, 'C4L': 1, 'C4S': 1},
    '<AL-BS-2>': {'C4A': 1, 'C4B':1, 'C4L': 1, 'C4S': 1},
    '<AL-BS-3>': {'C4A': 1, 'C4B':1, 'C4L': 1, 'C4S': 1},
    '<AL-BS-4>': {'C4A': 1, 'C4B':1, 'C4L': 1, 'C4S': 1},
    '<AL-BS-5>': {'C4A': 1, 'C4B':1, 'C4L': 1, 'C4S': 1},
    '<AL-BS-other>': {'C4A': 1, 'C4B':1, 'C4L': 1, 'C4S': 1},
    '<BS>': {'C4A': 0, 'C4B':1, 'C4L': 0, 'C4S': 1},
    '<AL-AS-BL>': {'C4A': 2, 'C4B':1, 'C4L': 2, 'C4S': 1},
    '<AL-BL-1>': {'C4A': 1, 'C4B': 1, 'C4L': 2, 'C4S': 0},
    '<AL-BL-2>': {'C4A': 1, 'C4B': 1, 'C4L': 2, 'C4S': 0},
    '<AL-BL-3>': {'C4A': 1, 'C4B': 1, 'C4L': 2, 'C4S': 0},
    '<AL-BL-other>': {'C4A': 1, 'C4B': 1, 'C4L': 2, 'C4S': 0},
    '<AL-AL-BL>':{'C4A': 2, 'C4B': 1, 'C4L': 3, 'C4S': 0},
    '<BL>': {'C4A': 0, 'C4B':1, 'C4L': 1, 'C4S': 0},
    '<AL-AS-BS-BS>': {'C4A': 2, 'C4B':2, 'C4L': 1, 'C4S': 3},
    '<AL>': {'C4A': 1, 'C4B':0, 'C4L': 1, 'C4S': 0},
    '<AL-BL-BL>': {'C4A': 1, 'C4B':2, 'C4L': 3, 'C4S': 0},
    '<AL-BS-BS>': {'C4A': 1, 'C4B':2, 'C4L': 1, 'C4S': 2},
    '<AL-AL-BS>': {'C4A': 2, 'C4B':1, 'C4L': 2, 'C4S': 1},
    '<AL-AL-AL>': {'C4A': 3, 'C4B':0, 'C4L': 3, 'C4S': 0},
    '<AL-BL-BS>': {'C4A': 1, 'C4B':2, 'C4L': 2, 'C4S': 1},
    '<BL-BS>': {'C4A': 0, 'C4B':2, 'C4L': 1, 'C4S': 1}
}


# Initialize dictionaries for counts and dosages
c4_isotype_counts = {sample: {'C4A': 0, 'C4B': 0, 'C4L': 0, 'C4S': 0} for sample in dosages_info.columns[3:]}
c4_isotype_dosages = {sample: {'C4A': 0.0, 'C4B': 0.0, 'C4L': 0.0, 'C4S': 0.0} for sample in dosages_info.columns[3:]}

# Process dosage values and counts
for _, row in dosages_info.iterrows():
    haplotype = row['ALT']
    dr2 = row['DR2']
    if haplotype in haplotypes and dr2 >= 0.3:
        for sample in dosages_info.columns[3:]:
            dosage = row[sample]
            if dosage != 0 and not pd.isnull(dosage):
                for isotype, count in haplotypes[haplotype].items():
                    c4_isotype_counts[sample][isotype] += count
                    c4_isotype_dosages[sample][isotype] += dosage * count

# Convert to DataFrame and save outputs
c4_isotype_counts_df = pd.DataFrame(c4_isotype_counts).transpose().reset_index()
c4_isotype_counts_df.rename(columns={'index': 'Sample_ID'}, inplace=True)

c4_isotype_dosages_df = pd.DataFrame(c4_isotype_dosages).transpose().reset_index()
c4_isotype_dosages_df.rename(columns={'index': 'Sample_ID'}, inplace=True)

c4_isotype_counts_df.to_csv(counts_file, index=False, sep='\t')
c4_isotype_dosages_df.to_csv(dosages_file, index=False, sep='\t')

print(f"\nProcessing complete. Output files:")
print(f"\nC4 isotype counts: {counts_file}"  )
print(f"\nC4 isotype dosages: {dosages_file}" )

#sum_dosages = c4_isotype_dosages_df_t.sum(axis=1)
#print("Sum of dosages for each sample:")
#print(sum_dosages)







