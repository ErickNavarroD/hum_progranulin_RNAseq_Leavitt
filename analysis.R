# Loadpackages ----
library(tidyverse)
library(DESeq2)

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

## read data ----

### metadata ----
metadata <- read.table(
  "data/metadata.tsv", header = FALSE,
  col.names = c("sample_id")
) %>% 
  mutate(
    condition = str_remove_all(sample_id, '\\d+'),
    condition = factor(
      condition,
      levels = c(
        'WT', 'Het', 'GhKO'
      )
    )
  ) %>% 
  column_to_rownames("sample_id")

### counts data ----
counts_data <- read.table(
  counts_data_path, header = TRUE, sep = "\t"
) %>% 
  dplyr::select("gene_id", all_of(rownames(metadata))) %>% 
  column_to_rownames("gene_id")

### make sure they match ----
stopifnot(
  all.equal(colnames(counts_data), rownames(metadata))
)

# DE analysis ----

## fit models ----
make_dir("tmp")
if (!file.exists("tmp/deseq2.rds")) {
  dds <- DESeqDataSetFromMatrix(countData = counts_data,
                                colData = metadata,
                                design = ~ condition) %>% 
    DESeq(betaPrior = TRUE)
  saveRDS(dds, "tmp/deseq2.rds")
} else {
  dds <- readRDS("tmp/deseq2.rds")
}




## retrieve results ----

comparisons <- list(
  het_vs_wt = c("condition", "Het", "WT"),  # log2FC is log2(Het/WT)
  ghko_vs_wt = c("condition", "GhKO", "WT"),
  ghko_vs_het = c("condition", "GhKO", "Het")
)

all_results <- imap_dfr(
  comparisons, 
  function(.comparison, .comparison_name) {
    
    res <- results(dds, contrast = .comparison, tidy = TRUE) %>% 
      dplyr::rename(gene_id := row) %>% 
      dplyr::filter(!is.na(padj)) %>% 
      mutate(comparison = .comparison_name) %>% 
      arrange(padj, -abs(log2FoldChange), -baseMean)
    
    
    output_dir <- str_glue("output/diff_expression_tables/{.comparison_name}")
    make_dir(output_dir)
    
    # write all results for this comparison 
    # (except for genes with NA in padj column)
    write_tsv(
      res,
      str_glue("{output_dir}/{.comparison_name}_all_genes.tsv")
    )
    
    # write significant calls
    write_tsv(
      res %>% filter(padj < 0.05),
      str_glue("{output_dir}/{.comparison_name}_significant_genes.tsv")
    )
    
    # write top ~5% calls
    # might be greater than 5% if there are ties in padj
    write_tsv(
      res %>% 
        filter(padj < 0.05) %>% 
        slice_min(padj, prop = .05),
      str_glue("{output_dir}/{.comparison_name}_top5_genes.tsv")
    )
    
    return(res)
  }
  
)

## look at Grn gene ----

data.frame(
  y = counts(dds, normalized=T)["Grn",],
  x = dds$condition
) %>% 
  ggplot(aes(x, y)) +
  ggbeeswarm::geom_quasirandom(width = .2) +
  ylim(0, 1500) +
  labs(x = NULL, y = "Normalized Counts",
       title = "Grn")
