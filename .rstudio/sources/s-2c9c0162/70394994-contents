# Load helper functions ----
source("code/utils.R")

# Merge counts data (if needed) ----

counts_data_path <- "data/counts.tsv"

if (!file.exists(counts_data_path)) {
  merge_counts_data(
    filepath_file = "data/filepaths.txt",
    output_path = counts_data_path
  )
  
}

counts_data_path <- read.table(
  counts_data_path, header = TRUE, sep = "\t"
)

# Differential expression analysis ----
