merge_counts_data <- function(
  filepath_file,
  output_path
) {
  filepath_vector <- readLines(filepath_file)
  df_list <-lapply(filepath_vector, function(filepath) {
    .df <- read.table(filepath)
    colnames(.df)[1] <- "gene_id"
    .sample <- strsplit(filepath, "/")[[1]][8]
    .sample <- strsplit(.sample, "_")[[1]][1]
    colnames(.df)[2] <- .sample
    return(.df)
  })
  
  counts_data <- Reduce(
    function(x,y) merge(x = x, y = y, by = "gene_id"), 
    df_list
  )
  
  write.table(
    counts_data,
    output_path,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
}

make_dir <- function(.path) {
  if (!dir.exists(.path)) {
    dir.create(.path, recursive = T)
  }
}
