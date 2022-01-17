# Leavitt collab - 2nd batch

Full analysis pipeline is in `analysis.R` file. 

## Directories

* `data/` all data necessary to replicate the analysis.
* `code/` auxiliary R functions.
* `output/` all output intended for final researcher's use.
  * `output/diff_expression_tables/` tables with results from differential abundance analysis.
  * `output/figures/` generated figures.
* `tmp/` intermediary objects intended solely to ease computing. 

## Files

* `data/counts.tsv` counts from merged Illumina/Basespace files. If absent, the first step in `analysis.R` creates this file, which should be run inside Sockeye.
* `data/filepaths.txt` paths for counts files for individual samples in Sockeye disk. The first step in `analysis.R` consults these paths to generate the above-cited `counts.tsv` file.

## Methods

All statistical analyses were carried out using R (v. 4.1.2) (R Core Team, 2021). Expression counts
were retrieved directly from BaseSpace application (Illumina, 2021). Distributions of read counts from each
sample were visually inspected before and after normalization to detect any sample with
uncommon overall expression profile. We performed differential expression (DE) analysis
using the DESeq2 R package (v. 1.34.0) (Love et al., 2014). P-values were adjusted using 
the Benjamini-Hochberg procedure to control False-discovery Rates at 5% (Benjamini & Hochberg, 1995). Gene set and enrichment analyses were carried out using the clusterProfiler R package (v. 4.2.2) (Wu et al., 2021). We used the tidyverse R package suit (v. 1.3.1) (Wickham et al., 2019) for general data wrangling and visualization.


## References

* R Core Team (2021). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.

* Love, M.I., Huber, W., Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550. https://doi.org/10.1186/s13059-014-0550-8

* Benjamini, Y., & Hochberg, Y. (1995). Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. In Journal of the Royal Statistical Society: Series B (Methodological) (Vol. 57, Issue 1, pp. 289–300). Wiley. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x

* Wu T, Hu E, Xu S, Chen M, Guo P, Dai Z, Feng T, Zhou L, Tang W, Zhan L, Fu x, Liu S, Bo X, Yu G. (2021). “clusterProfiler 4.0: A universal enrichment tool for interpreting omics data.” The Innovation, 2(3), 100141. https://doi.org/10.1016/j.xinn.2021.100141

* Wickham, H., Averick, M., Bryan, J., Chang, W., McGowan, L., François, R., Grolemund, G., Hayes, A., Henry, L., Hester, J., Kuhn, M., Pedersen, T., Miller, E., Bache, S., Müller, K., Ooms, J., Robinson, D., Seidel, D., Spinu, V., … Yutani, H. (2019). Welcome to the Tidyverse. In Journal of Open Source Software (Vol. 4, Issue 43, p. 1686). The Open Journal. https://doi.org/10.21105/joss.01686
