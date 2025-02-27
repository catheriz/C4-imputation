# Install and load necessary packages
library(tidyverse)


data <- read.table('dosages.txt',sep='\t',head=T,check.names = FALSE)

process_info <- function(info) {
  strsplit(info, ",")
}

# List of columns to process
columns_to_process <- names(data)

# Process each column to split comma-separated values and unnest
for (col in columns_to_process) {
    data <- data %>%
      mutate(!!col := strsplit(as.character(!!sym(col)), ","))
}

# Unnest all columns
processed_data <- data %>%
  unnest(cols = all_of(columns_to_process))

# View the processed data
#print(processed_data)

# Save the processed data
write.table(processed_data, "processed_data.txt", row.names = F,col.names=T,sep='\t',quote=F)
