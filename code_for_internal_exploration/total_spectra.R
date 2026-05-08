#!/usr/bin/env Rscript
# For each indel classification (83, 89, 476), sum the per-tumor spectra
# (row sums) into a single "Total" spectrum and plot all three on one
# page in Manuscript_data/finalized_cap9/total_spectra.pdf.

library(mSigPlot)
library(ggplot2)
library(gridExtra)

types <- c("83", "89", "476")

plot_fns <- list(
  "83"  = mSigPlot::plot_ID83,
  "89"  = mSigPlot::plot_ID89,
  "476" = mSigPlot::plot_ID476
)

plot_args <- list(
  "83"  = list(block_label_cex = 0.7, axis_text_x_cex = 0.5, axis_title_y_cex = 0.1),
  "89"  = list(block_label_cex = 1.5, axis_text_x_cex = 0.4, axis_title_y_cex = 0.1),
  "476" = list()
)

plots <- lapply(types, function(type) {
  spec_file <- here::here(
    "Manuscript_data/finalized_cap9",
    paste0("liu_et_al_", type, "_spectra.tsv")
  )
  spec <- read.delim(spec_file, sep = "\t", row.names = 1, check.names = FALSE)
  total <- data.frame(Total = rowSums(spec), row.names = rownames(spec),
                      check.names = FALSE)
  do.call(plot_fns[[type]], c(list(total), plot_args[[type]]))
})

out_path <- here::here("Manuscript_data/finalized_cap9", "total_spectra.pdf")
pdf(out_path, width = 8.5, height = 11)
grid.arrange(grobs = plots, ncol = 1, nrow = 3)
dev.off()
message("Wrote ", out_path)
