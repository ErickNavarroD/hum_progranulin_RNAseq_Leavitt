# Loadpackages ----
library(tidyverse)
library(DESeq2)
library(here)
library(ggrepel)
library(EnhancedVolcano)
library(pheatmap)
library(ggvenn)
#if (!require(devtools)) install.packages("devtools")
#devtools::install_github("yanlinlin82/ggvenn")
library(clusterProfiler)
library(org.Mm.eg.db)
# Load helper functions ----
source("code/utils.R")

# Load data ----

counts_data_path <- here("data","counts.tsv")

## merge counts files if needed ----
if (!file.exists(counts_data_path)) {
  merge_counts_data(
    filepath_file = "data/filepaths.txt",
    output_path = counts_data_path
  )
  
}

## read data ----

### metadata ----
#Create the metadata file
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
  column_to_rownames("sample_id") %>% 
  mutate(genotype = factor(c(rep("Het",8),rep("Hom",5)), # Condition that will be used for the mixed-samples analysis (mouse&human vs wt)
                           levels = c("Hom", "Het")))

### counts data ----
counts_data <- read.table(counts_data_path, header = TRUE, sep = "\t") %>% 
  dplyr::select("gene_id", all_of(rownames(metadata))) %>% #Sort the columns
  column_to_rownames("gene_id")

### make sure metadata and counts table match ----
stopifnot(
  all.equal(colnames(counts_data), rownames(metadata))
)

# DE analysis ----

## fit models ----
# Create the dds object or load it if already existing
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

#Make a list with the comparisons we want to make
comparisons <- list(
  het_vs_wt = c("condition", "Het", "WT"),  # log2FC is log2(Het/WT)
  ghko_vs_wt = c("condition", "GhKO", "WT"),
  ghko_vs_het = c("condition", "GhKO", "Het")
)

#Create an object with the results for each comparison
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

# Plotting ----

## PCA
make_dir(here("output", "figures"))

vsd = vst(object = dds, blind = TRUE)
jpeg(file = here("output", "figures", "PCA.jpeg"))
plotPCA(vsd, intgroup = c("condition")) + 
  theme_classic()+
  scale_colour_manual(values= c("#999999", "#E69F00", "#56B4E9"))+
  geom_text_repel(aes(label = row.names(metadata)))
dev.off()

## Volcano plots

for (contrast in names(comparisons)){
  jpeg(file = here("output","figures",paste(contrast,"_VolcPlot.jpeg",sep = "")))
  print(all_results %>% 
    filter(comparison == contrast) %>% 
    EnhancedVolcano(lab = .[,1],
                    x = "log2FoldChange", 
                    FCcutoff = 0, 
                    y = "padj",
                    title = contrast,
                    subtitleLabSize = 0.1,
                    legendLabels=c('NS', 
                                   "log2FC",
                                   'p adj',
                                   "p adj & log2FC"))+
    ylab(expression(-log[10]~p~adj))) 
    dev.off()
}

  


## look at the expression of the Grn gene ----

jpeg(file = here("output","figures","Grn_expression.jpeg"))
data.frame(
  y = counts(dds, normalized=T)["Grn",],
  x = dds$condition
) %>% 
  ggplot(aes(x, y)) +
  ggbeeswarm::geom_quasirandom(width = .2) +
  ylim(0, 1500) +
  labs(x = NULL, y = "Normalized Counts",
       title = "Grn") +
  theme_classic()+
  ggtitle("Grn expression")
dev.off()

## Heatmaps ----

#Get normalized counts that will be used to plot the heatmap
norm_counts =
  vsd %>% 
  assay() %>% as.data.frame() %>% 
  rownames_to_column() %>% 
  filter(rowname %in% ( #Select the genes that were DE in any comparison
    all_results %>% 
      filter(padj < 0.05) %>% 
      pull(gene_id) %>% 
      unique())) %>% column_to_rownames(var = "rowname")

#Plot the heatmap
jpeg(file = here("output", "figures", "heatmap_unclustered.jpeg"))
pheatmap(norm_counts, scale = "row",
         cluster_cols = F,
         main = "DE genes in any comparison (1382)", #Poner el titulo
         show_rownames = F,
         show_colnames = T, 
         treeheight_col = 20,
         treeheight_row = 0, #Remove row dendogram
         annotation_col = metadata %>% select(condition),
         annotation_colors = list(condition = c(WT = "#999999", Het = "#E69F00", GhKO = "#56B4E9")))
dev.off() 

## Plot the Venn diagram ----

#Plot the Number of DE genes
jpeg(file = here("output", "figures","DE_genes_numbers.jpeg"))
all_results %>% 
  filter(padj < 0.05) %>% 
  mutate(Change = case_when(log2FoldChange > 0 ~ "Overexpressed" ,
                            log2FoldChange < 0 ~ "Underexpressed")) %>% 
  ggplot(aes(x = comparison, fill = Change)) +
  geom_bar()+
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = 0)+
  theme_classic()+
  ylab("Number of genes")+
  ggtitle("DE genes per comparison")+
  scale_x_discrete("") #Remove x axis title
dev.off()  
  

#Plot the Venn diagram
venn_info <- list(
  Het = all_results %>% 
    filter(padj < 0.05,
           comparison == "het_vs_wt") %>% 
    pull(gene_id),
  GhKO = all_results %>% 
    filter(padj < 0.05,
           comparison == "ghko_vs_wt") %>% 
    pull(gene_id)
)

jpeg(file = here("output", "figures", "Venn_Het_GhKO.jpeg"))
ggvenn(venn_info,stroke_size = 0.5)
dev.off()

## Enrichment analysis ----
universo = unique(all_results$gene_id)
number_enriched_terms = tibble()
for (contrast in names(comparisons)){
  
  subset_comp = all_results %>% 
    filter(comparison ==contrast,
           padj < 0.05)
  
  enrich_GO_res = enrichGO(gene = as.character(na.omit(subset_comp$gene_id)),
                     keyType = "SYMBOL",
                     OrgDb = 'org.Mm.eg.db',
                     ont="BP", pvalueCutoff=0.05,qvalueCutoff = 0.5,
                     universe = as.character(na.omit(universo)),
                     readable = F)
  
  jpeg(file = here("output","figures",paste(contrast,"_GObp_dotplot.jpeg",sep = "")))
  print(
    clusterProfiler::dotplot(enrich_GO_res, showCategory = 25, title = contrast )
  ) 
  dev.off()
  
  enriched_terms = enrich_GO_res@result %>% 
    filter(p.adjust < 0.05) %>% 
    nrow()
  
  number_enriched_terms = c(number_enriched_terms, enriched_terms )
}

jpeg(file = here("output", "figures", "Number_enriched_terms.jpeg"))
number_enriched_terms %>% 
  as.tibble() %>% 
  rename(enriched_terms = value) %>% 
  mutate(contrast = names(comparisons)) %>% 
  ggplot(aes(x = contrast, y = enriched_terms))+
  geom_col(aes(fill = contrast), alpha = 0.5)+
  theme_classic()+
  ylab("Enriched terms")+
  guides(fill = "none") +
  scale_x_discrete("")
dev.off()  

## Het + GhKO vs WT ----
dds_genotype <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = metadata,
                              design = ~ genotype) %>% 
  DESeq(betaPrior = TRUE)

res_genotype = results(dds_genotype, tidy = TRUE) %>% 
  dplyr::rename(gene_id := row) %>% 
  dplyr::filter(!is.na(padj)) %>% 
  arrange(padj, -abs(log2FoldChange), -baseMean)

res_gen_signif = res %>% 
  filter(padj < 0.05)

## Write the results ----
output_dir <- str_glue("output/diff_expression_tables/ghkoHet_vs_wt")
make_dir(output_dir)
# write all results for this comparison 
# (except for genes with NA in padj column)
write_tsv(
  res_genotype,
  str_glue("{output_dir}/ghkoHet_vs_wt_all_genes.tsv")
)

# write significant calls
write_tsv(
  res_genotype %>% filter(padj < 0.05),
  str_glue("{output_dir}/ghkoHet_vs_wt_significant_genes.tsv")
)

# write top ~5% calls
# might be greater than 5% if there are ties in padj
write_tsv(
  res_genotype %>% 
    filter(padj < 0.05) %>% 
    slice_min(padj, prop = .05),
  str_glue("{output_dir}/ghkoHet_vs_wt_top5_genes.tsv")
)

## Volcano plot ----
jpeg(file = here("output","figures","mixed_genotypes_volcanoplot.jpeg"))
res_genotype %>% 
  EnhancedVolcano(lab = .[,1],
                  x = "log2FoldChange", 
                  FCcutoff = 0, 
                  y = "padj",
                  title = "GhKO + Het vs WT",
                  subtitleLabSize = 0.1,
                  legendLabels=c('NS', 
                                 "log2FC",
                                 'p adj',
                                 "p adj & log2FC"))+
  ylab(expression(-log[10]~p~adj)) 
dev.off()

#Barplot of number of DE genes
jpeg(file = here("output","figures","mixed_genotypes_numberDEgenes.jpeg"))
res_genotype %>% 
  filter(padj < 0.05) %>% 
  mutate(Change = case_when(log2FoldChange > 0 ~ "Overexpressed" ,
                            log2FoldChange < 0 ~ "Underexpressed"),
         comparison = "het+GhKO_vs_WT") %>% 
  ggplot(aes(x = comparison, fill = Change)) +
  geom_bar()+
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = 0)+
  theme_classic()+
  ylab("Number of genes")+
  ggtitle("DE genes")+
  scale_x_discrete("") #Remove x axis title
dev.off()

#Comparison with a Venn diagram with previous results
jpeg(file = here("output","figures","mixed_genotypes_Venn.jpeg"))
venn_info <- list(
  Het = all_results %>% 
    filter(padj < 0.05,
           comparison == "het_vs_wt") %>% 
    pull(gene_id),
  GhKO = all_results %>% 
    filter(padj < 0.05,
           comparison == "ghko_vs_wt") %>% 
    pull(gene_id),
  GhKO_Het = res_genotype %>% 
    filter(padj < 0.05) %>% 
    pull(gene_id)
)
ggvenn(venn_info,stroke_size = 0.5)
dev.off()

## GO enrichment analysis
universo_genotype = unique(res_genotype$gene_id)

enrich_GO_res_genotype = enrichGO(gene = as.character(na.omit(res_gen_signif$gene_id)),
                         keyType = "SYMBOL",
                         OrgDb = 'org.Mm.eg.db',
                         ont="BP", pvalueCutoff=0.05,qvalueCutoff = 0.5,
                         universe = as.character(na.omit(universo_genotype)),
                         readable = F)

jpeg(file = here("output","figures","mixed_genotypes_GOenrich.jpeg"))
clusterProfiler::dotplot(enrich_GO_res, showCategory = 25, title = contrast ) 
dev.off()
#No enriched terms were found. 
