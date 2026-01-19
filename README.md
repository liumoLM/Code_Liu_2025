# Code and Data for Liu et al.

File Manuscript_data/Koh_signatures.tsv came from Supplementary Table 10 in Koh et al., 2024

To render single-file html files for each signature, first render vignette.qmd, then at bash
command line run

`Rscript render_separate_pages.R`

library(quarto)
setwd("vignette")
quarto_render("vignette.qmd", output_format = "html")
