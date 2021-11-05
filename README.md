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
