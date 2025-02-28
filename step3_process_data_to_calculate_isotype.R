# Install and load necessary packages
library(tidyverse)

# Function to display usage
usage <- function() {
  cat("Usage: Rscript step3_process_data_to_calculate_isotype.R<input_file> <output_file>\n")
  cat("Example: Rscript step3_process_data_to_calculate_isotype.R dosages.txt processed_data.txt\n")
  quit(status = 1)
}

# Check if correct number of arguments are provided
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  usage()
}

# Read input and output file names from command-line arguments
input_file <- args[1]
output_file <- args[2]

# Read the data
data <- read.table(input_file, sep = '\t', header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

# Vectorized: Apply strsplit to all columns at once
processed_data <- data %>%
  mutate(across(everything(), ~ strsplit(as.character(.), ","))) %>%
  unnest(cols = everything())

# Save the processed data
write.table(processed_data, output_file, row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)

cat("Processing complete. Output saved to:", output_file, "\n")
