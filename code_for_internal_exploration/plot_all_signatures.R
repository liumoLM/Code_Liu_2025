#!/usr/bin/env Rscript
# Plot all finalized signatures of a given ID type to a multi-page PDF.
#
# Layout (8.5 x 11 inch pages, portrait):
#   83-type and 89-type:  2 columns x 4 rows per page
#   476-type:             1 column  x 4 rows per page
#
# Arguments:
#   --type  Signature type to plot: "83", "89", or "476" (required)
#
# Output:
#   plot_output/all_{type}_signatures.pdf
#
# Usage:
#   Rscript plot_all_signatures.R --type 83
#   Rscript plot_all_signatures.R --type 89
#   Rscript plot_all_signatures.R --type 476

library(argparser)
library(mSigPlot)
library(ggplot2)
library(gridExtra)

p <- arg_parser("Plot all signatures to PDF")
p <- add_argument(
  p,
  "--type",
  help = "Signature type: 83, 89, or 476",
  type = "character"
)
args <- parse_args(p)

plot_fn <- switch(
  args$type,
  "83" = mSigPlot::plot_ID83,
  "89" = mSigPlot::plot_ID89,
  "476" = mSigPlot::plot_ID476,
  stop("--type must be 83, 89, or 476")
)

sig_file <- here::here(
  "Manuscript_data/finalized_cap9",
  paste0("liu_et_al_", args$type, "_signatures.tsv")
)

sigs <- read.delim(sig_file, sep = "\t", row.names = 1, check.names = FALSE)

n_sigs <- ncol(sigs)

if (args$type == "476") {
  ncols <- 1
  nrows <- 4
} else if (args$type == "89") {
  ncols <- 2
  nrows <- 5
} else {
  ncols <- 2
  nrows <- 4
}
per_page <- ncols * nrows

out_path <- here::here(
  "plot_output",
  paste0("all_", args$type, "_signatures.pdf")
)
# Extra plot arguments per type
plot_args <- if (args$type == "83") {
  list(block_label_cex = 0.7, axis_text_x_cex = 0.5, axis_title_y_cex = 0.1)
} else if (args$type == "89") {
  list(block_label_cex = 1.5, axis_text_x_cex = 0.4, axis_title_y_cex = 0.1)
} else {
  list()
}

pdf(out_path, width = 8.5, height = 11)

for (start in seq(1, n_sigs, by = per_page)) {
  end <- min(start + per_page - 1, n_sigs)
  plots <- lapply(start:end, function(i) {
    do.call(plot_fn, c(list(sigs[, i, drop = FALSE]), plot_args))
  })
  grid.arrange(grobs = plots, ncol = ncols, nrow = nrows)
}

dev.off()
message("Wrote ", out_path)
