# Load helper functions ----
source("code/utils.R")

# Load data ----

counts_data_path <- "data/counts.tsv"

## merge counts files if needed ----
if (!file.exists(counts_data_path)) {
  merge_counts_data(
    filepath_file = "data/filepaths.txt",
    output_path = counts_data_path
  )
  
}

## read counts and metadata ----
counts_data_path <- read.table(
  counts_data_path, header = TRUE, sep = "\t"
)

metadata <- read.table(
  "data/metadata.tsv", header = FALSE,
  col.names = c("sample_id")
)

# DE analysis ----
