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

plot_gsea <- function(df, N = 25) {
  df %>% 
    arrange(p.adjust) %>% 
    head(N) %>% 
    mutate(
      Description = fct_reorder(
        Description,
        enrichmentScore
      ),
      x = map_dbl(core_enrichment, ~ {
        str_split(.x, fixed("/")) %>% 
          (
            function(k) length(k[[1]])
          )
      }),
      .gene_ratio = x/setSize
    ) %>% 
    ggplot(aes(enrichmentScore, Description)) +
    geom_point(
      aes(color = p.adjust, size = .gene_ratio)
    ) +
    theme_bw() +
    labs(x = "Enrichment Score",
         y = NULL, 
         color = "p adjusted",
         size = "Gene Ratio") +
    scale_color_viridis_b(
      breaks = scales::pretty_breaks(7)
    ) +
    theme(legend.key.height = unit(1.2, "cm"))
}

