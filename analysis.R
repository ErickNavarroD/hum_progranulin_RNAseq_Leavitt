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
all_results_file <- here("output/diff_expression_tables/all_results.tsv")
if (!file.exists(all_results_file)) {
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
  
  write_tsv(all_results, all_results_file)
  
} else {
  
  all_results <- read_tsv(all_results_file)
  
}
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


## Differentially expressed genes ----

heatmap_columns_ordered <- c(
  paste0("WT", 1:5),
  paste0("Het", 1:4),
  paste0("GhKO", 1:4)
)

#Get normalized counts that will be used to plot the heatmap
norm_counts =
  vsd %>% 
  assay() %>% as.data.frame() %>% 
  rownames_to_column() %>% 
  filter(rowname %in% ( #Select the genes that were DE in any comparison
    all_results %>% 
      filter(padj < 0.05) %>% 
      pull(gene_id) %>% 
      unique())) %>% column_to_rownames(var = "rowname") %>% 
  dplyr::select(heatmap_columns_ordered)

#Plot the heatmap
jpeg(file = here("output", "figures", "heatmap_unclustered.jpeg"),
     res = 300, width = 1600, height = 1600)
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

## Selected genes ----
genes_of_interest <- c( 
  "Apoa1", "Apoa2", "Apoa5", "Apob", "Apoc1", "Apoc3", 
  "Aloxe3", "Osbpl5", "Cyp2d9", "Serpina1b", "Serpina1c", "Serpina1d"
)
norm_counts_selected = norm_counts[genes_of_interest, ]


#Plot the heatmap
jpeg(file = here("output", "figures", "heatmap_selected_genes.jpeg"),
     res = 300, width = 1600, height = 1600)
pheatmap(norm_counts_selected,
         scale = "row",
         cluster_cols = F,
         main = "Selected DE genes", 
         show_rownames = T,
         show_colnames = T, 
         treeheight_col = 20,
         treeheight_row = 0, #Remove row dendogram
         annotation_col = metadata %>% select(condition),
         annotation_colors = list(condition = c(WT = "#999999", Het = "#E69F00", GhKO = "#56B4E9")))
dev.off() 

## Enrichment analysis ----
universo = unique(all_results$gene_id)

for (fc_cutoff in c(0,.5,1)){ #do the analysis for each subontology
  dir.create(here("output",
                  "figures",
                  str_glue("exploratory_FCcuts_",as.character(fc_cutoff)))
  )
  ## Plot the Venn diagram ----
  #Plot the Number of DE genes
  jpeg(file = here("output",
                   "figures",
                   str_glue("exploratory_FCcuts_",as.character(fc_cutoff)),
                   str_glue("fc_",as.character(fc_cutoff),"_","DE_genes_numbers.jpeg")))
  print(
    all_results %>% 
      filter(padj < 0.05,
             abs(log2FoldChange) > fc_cutoff) %>% 
      mutate(Change = case_when(log2FoldChange > 0 ~ "Overexpressed" ,
                                log2FoldChange < 0 ~ "Underexpressed")) %>% 
      ggplot(aes(x = comparison, fill = Change)) +
      geom_bar()+
      geom_text(stat = "count", aes(label = after_stat(count)), vjust = 0)+
      theme_classic()+
      ylab("Number of genes")+
      ggtitle("DE genes per comparison")+
      scale_x_discrete("")
  ) #Remove x axis title
  dev.off()  
  
  
  #Plot the Venn diagram
  venn_info <- list(
    Het = all_results %>% 
      filter(padj < 0.05,
             abs(log2FoldChange) > fc_cutoff,
             comparison == "het_vs_wt") %>% 
      pull(gene_id),
    GhKO = all_results %>% 
      filter(padj < 0.05,
             abs(log2FoldChange) > fc_cutoff,
             comparison == "ghko_vs_wt") %>% 
      pull(gene_id)
  )
  
  jpeg(file =  here("output",
                    "figures",
                    str_glue("exploratory_FCcuts_",as.character(fc_cutoff)),
                    str_glue("fc_",as.character(fc_cutoff),"_","Venn_Het_GhKO.jpeg")))
  print(ggvenn(venn_info,stroke_size = 0.5))
  dev.off()
  #
  for (subontology in c("BP","MF","CC")){ #Try an enrichment analysis for each  subontology
    number_enriched_terms = c()
    for (contrast in names(comparisons)){ #Do it for each comparison
      
      subset_comp = all_results %>% 
        filter(comparison ==contrast,
               padj < 0.05,
               abs(log2FoldChange) > fc_cutoff)
      
      enrich_GO_res = enrichGO(gene = as.character(na.omit(subset_comp$gene_id)),
                               keyType = "SYMBOL",
                               OrgDb = 'org.Mm.eg.db',
                               ont=subontology, pvalueCutoff=0.05,qvalueCutoff = 0.5,
                               universe = as.character(na.omit(universo)),
                               readable = F)
      
      
      #Append the number of enriched terms to the object to plot later on the number per comparison in a barplot
      if (!is.null(enrich_GO_res)){
        enriched_terms = enrich_GO_res@result %>% 
          filter(p.adjust < 0.05) %>% 
          nrow()} else{
            enriched_terms = 0
          }
      
      number_enriched_terms = c(number_enriched_terms, enriched_terms )
      
      #Plot the enriched terms if there is any
      if (enriched_terms > 0){
        jpeg(file = here("output",
                         "figures",
                         str_glue("exploratory_FCcuts_",as.character(fc_cutoff)),
                         str_glue(contrast,"_fc_",as.character(fc_cutoff),"_",subontology,"_dotplot.jpeg",sep = "")))
        print(clusterProfiler::dotplot(enrich_GO_res, showCategory = 25, title = contrast ))
        dev.off()
        
      }
      
    }
    
    jpeg(file = here("output",
                     "figures",
                     str_glue("exploratory_FCcuts_",as.character(fc_cutoff)),
                     str_glue(subontology,"_fc_",as.character(fc_cutoff),"_","Number_enriched_terms.jpeg")))
    
    print(number_enriched_terms %>% 
            as_tibble() %>% 
            rename(enriched_terms = value) %>% 
            mutate(contrast = names(comparisons)) %>% 
            ggplot(aes(x = contrast, y = enriched_terms))+
            geom_col(aes(fill = contrast), alpha = 0.5)+
            theme_classic()+
            ylab("Enriched terms")+
            guides(fill = "none") +
            scale_x_discrete(""))
    dev.off()
  }
}

# Enrichment of the 302 common dysregulated genes in the het vs wt and ghko vs wt comparisons    

venn_info <- list(
  Het = all_results %>% 
    filter(padj < 0.05,
           abs(log2FoldChange) > 0,
           comparison == "het_vs_wt") %>% 
    pull(gene_id),
  GhKO = all_results %>% 
    filter(padj < 0.05,
           abs(log2FoldChange) > 0,
           comparison == "ghko_vs_wt") %>% 
    pull(gene_id)
)

comm_dys_genes = intersect(venn_info$Het , venn_info$GhKO)

dir.create(here("output",
                "figures",
                "comm_dys_genes")
)

for (subontology in c("BP","MF","CC")){ #Try an enrichment analysis for each  subontology
  number_enriched_terms = c()
  
  enrich_GO_res = enrichGO(gene = comm_dys_genes,
                           keyType = "SYMBOL",
                           OrgDb = 'org.Mm.eg.db',
                           ont=subontology, pvalueCutoff=0.05,qvalueCutoff = 0.5,
                           universe = as.character(na.omit(universo)),
                           readable = F)
  
  #Append the number of enriched terms to the object to plot later on the number per comparison in a barplot
  if (!is.null(enrich_GO_res)){
    enriched_terms = enrich_GO_res@result %>% 
      filter(p.adjust < 0.05) %>% 
      nrow()} else{
        enriched_terms = 0
      }
  
  number_enriched_terms = c(number_enriched_terms, enriched_terms )
  
  #Plot the enriched terms if there is any
  if (enriched_terms > 0){
    jpeg(file = here("output",
                     "figures",
                     "comm_dys_genes",
                     str_glue(subontology,"_dotplot.jpeg",sep = "")))
    print(clusterProfiler::dotplot(enrich_GO_res, showCategory = 25, title = contrast ))
    dev.off()
    
  }
  
}

jpeg(file = here("output",
                 "figures",
                 "comm_dys_genes",
                 str_glue("Number_enriched_terms.jpeg")))

print(number_enriched_terms %>% 
        as_tibble() %>% 
        rename(enriched_terms = value) %>% 
        mutate(contrast = names(comparisons)) %>% 
        ggplot(aes(x = contrast, y = enriched_terms))+
        geom_col(aes(fill = contrast), alpha = 0.5)+
        theme_classic()+
        ylab("Enriched terms")+
        guides(fill = "none") +
        scale_x_discrete(""))
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

# Check outliers ----

## normalized ----
jpeg("output/figures/sample_distributions_normalized.jpeg",
     res = 300, width = 1600, height = 1200)
counts(dds, normalized = TRUE) %>% 
  (\(x) log(x+1)) %>% 
  boxplot(las=2, ylab = "log[reads+1]", main = "Normalized")
dev.off()

## unnormalized ----
jpeg("output/figures/sample_distributions_unnormalized.jpeg",
     res = 300, width = 1600, height = 1200)
counts(dds, normalized = FALSE) %>% 
  (\(x) log(x+1)) %>% 
  boxplot(las=2, ylab = "log[reads+1]", main = "Not normalized")
dev.off()

# GSEA ----

## feature 1: numeric vector

for (contrast in names(comparisons)) {
  print(str_glue("Running GASEA for {contrast}"))
  res <- all_results %>% 
    filter(comparison == contrast) %>% 
    arrange(-log2FoldChange)
  geneList = res$log2FoldChange
  names(geneList) = as.character(res$gene_id)
  gsea_res <- gseGO(geneList = geneList,
                    OrgDb = org.Mm.eg.db,
                    ont = "ALL",
                    keyType = "SYMBOL",
                    minGSSize = 10,
                    eps = 1e-16,
                    maxGSSize = 500,
                    pvalueCutoff = 0.05,
                    verbose = FALSE)
  
  outdir <- str_glue("output/figures/gsea/{contrast}")
  make_dir(outdir)
  gsea_result_arranged <- gsea_res@result %>% 
    dplyr::filter(abs(enrichmentScore) > 0.5) %>% 
    arrange(p.adjust)
  write_tsv(gsea_result_arranged, 
            str_glue("{outdir}/{contrast}_gsea_stats.tsv"))
  
  n_sets <- 20
  .p <- dotplot(
    gsea_res, 
    showCategory = min(n_sets, nrow(gsea_result_arranged)),
    font.size = 9
  )
  ggsave(
    str_glue("{outdir}/{contrast}_gsea_dotplot_top{n_sets}_gene_ratio.jpeg"),
    .p, width = 9, height = 9
  )
  
  .p <- plot_gsea(gsea_result_arranged, N = n_sets)
  ggsave(
    str_glue("{outdir}/{contrast}_gsea_dotplot_top{n_sets}_enrich_score.jpeg"),
    .p, width = 12, height = 7
  )
  sets_to_plot <- gsea_result_arranged %>% 
    head(n_sets) %>% 
    pull(ID)
  
  top_sets_dir <- str_glue("{outdir}/top{n_sets}")
  make_dir(top_sets_dir)
  print(str_glue("Saving set plots"))
  for (set_id in sets_to_plot) {
    set_description <- gsea_result_arranged %>% 
      dplyr::filter(ID == set_id) %>% 
      pull(Description)
    .p <- gseaplot(gsea_res, 
                   geneSetID = set_id, 
                   title = str_to_title(set_description),
                   plot = FALSE) &
      theme(plot.title = element_text(size = 8))
    # .p <- .p & labs(
    #   title = paste0("\n", str_to_title(set_description))
    # )
    
    file_name <- paste0(
      str_remove(set_id, ":"), "_",
      str_replace_all(set_description, " |\\/", "-"),
      ".jpeg"
    )
    ggsave(
      str_glue("{top_sets_dir}/{file_name}"),
      .p, width = 8, height = 7.5
    )
  }
}


# Cluster comparisons ----

stopifnot(
  all(!near(all_results$log2FoldChange, 0))
)

all_results$direction <- ifelse(all_results$log2FoldChange > 0, "Up", "Down") 
all_results$comparison <- factor(
  as.character(all_results$comparison),
  levels = c('het_vs_wt', 'ghko_vs_wt', 'ghko_vs_het'),
  labels = c(
    latex2exp::TeX("$Grn^{+/-}\\ vs\\ Grn^{+/+}$"),
    latex2exp::TeX("$Grn^{-/-};GRN^{Tg}\\ vs \\ Grn^{+/+}$"),
    latex2exp::TeX("$Grn^{-/-};GRN^{Tg}\\ vs \\ Grn^{+/-}$")
  )
)

for (subontology in c("MF", "BP", "CC", "ALL")) {
  print(str_glue("Running enrichGO for subontology {subontology}"))  
  cluster_comp_mf = compareCluster(
    data = all_results %>% dplyr::filter(padj < .05),
    geneClusters = gene_id ~ direction + comparison, 
    fun = "enrichGO", 
    keyType = "SYMBOL",
    OrgDb = 'org.Mm.eg.db', 
    universe = unique(all_results$gene_id),
    ont = subontology
  )
  
  .plot <- dotplot(
    cluster_comp_mf,
    x = 'direction',
    showCategory = 5, 
    font.size = 9
  ) + 
    facet_wrap(
      ~comparison, 
      labeller = label_parsed
    ) +
    theme_bw(base_size = 18) +
    theme(legend.position = 'bottom',
          legend.key.width = unit(1, 'cm'),
          axis.text.y = element_text(size = 14)) +
    scale_y_discrete(
      labels = label_wrap_gen(width = 60)
    ) +
    scale_size_continuous(
      range = c(5, 12)
    )
  
  subontology <- str_to_lower(subontology)
  file_name <- str_glue(
    "output/stratified-go-analysis/stratified-go-enrichment-{subontology}.png"
  )
  ggsave(
    filename = file_name,
    plot = .plot,
    width = 16, height = 12
  )
}
