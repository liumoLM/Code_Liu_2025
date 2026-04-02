#!/usr/bin/env Rscript
# Plot all ID83 signatures to a multi-page PDF (2 columns x 5 rows per page)

library(mSigPlot)
library(ggplot2)
library(gridExtra)

sigs <- read.delim(
  here::here("Manuscript_data/finalized_cap9/liu_et_al_83_signatures.tsv"),
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

n_sigs <- ncol(sigs)
per_page <- 10  # 2 columns x 5 rows

pdf("/tmp/all_83_signatures.pdf", width = 20, height = 14)

for (start in seq(1, n_sigs, by = per_page)) {
  end <- min(start + per_page - 1, n_sigs)
  plots <- lapply(start:end, function(i) {
    mSigPlot::plot_ID83(sigs[, i, drop = FALSE])
  })
  grid.arrange(grobs = plots, ncol = 2, nrow = 5)
}

dev.off()
message("Wrote /tmp/all_83_signatures.pdf")
