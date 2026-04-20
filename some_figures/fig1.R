# fig1.R — 3-panel stacked indel-signature figure (ID83 / ID89 / ID476)
# Portrait letter PDF with editable text (cairo_pdf).

library(ggplot2)

# ---- Params (edit these to tune the figure) -------------------------------

sig_col_83  <- "C_ID23"
sig_col_89  <- "InsDel23"
sig_col_476 <- "InsDel23"

num_peaks   <- 4

base_size_83  <- 11
base_size_89  <- 11
base_size_476 <- 11

page_w      <- 8.5   # letter portrait width  (in)
page_h      <- 11    # letter portrait height (in)
plot_w      <- 6.5   # plot-block width       (in)
plot_h      <- 6.0   # plot-block height      (in) — ggrepel needs ~2 in per panel
margin_left <- 1     # from left edge         (in)
margin_top  <- 1     # from top edge          (in)

out_file <- here::here("some_figures", "fig1.pdf")

# ---- Read signatures ------------------------------------------------------

sig_dir <- here::here("Manuscript_data", "finalized_cap9")

sigs_83 <- read.delim(
  file.path(sig_dir, "liu_et_al_83_signatures.tsv"),
  row.names = 1, check.names = FALSE)

sigs_89 <- read.delim(
  file.path(sig_dir, "liu_et_al_89_signatures.tsv"),
  row.names = 1, check.names = FALSE)

sigs_476 <- read.delim(
  file.path(sig_dir, "liu_et_al_476_signatures.tsv"),
  row.names = 1, check.names = FALSE)

# ---- Build the three panels -----------------------------------------------

p83 <- mSigPlot::plot_ID83(
  sigs_83[, sig_col_83, drop = FALSE],
  num_peak_labels = num_peaks,
  base_size       = base_size_83)

p89 <- mSigPlot::plot_ID89(
  sigs_89[, sig_col_89, drop = FALSE],
  num_peak_labels = num_peaks,
  base_size       = base_size_89)

p476 <- mSigPlot::plot_ID476(
  sigs_476[, sig_col_476, drop = FALSE],
  num_peak_labels = num_peaks,
  base_size       = base_size_476)

# Strip outer plot margins so panels fill their 6.5" container horizontally
# and sit flush against each other vertically.
zero_margin <- ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0, "pt"))
p83  <- p83  + zero_margin
p89  <- p89  + zero_margin
p476 <- p476 + zero_margin

combined <- gridExtra::arrangeGrob(
  grobs = list(p83, p89, p476),
  ncol  = 1, nrow = 3)

# ---- Write PDF with plot block positioned in top printable area -----------

cairo_pdf(out_file, width = page_w, height = page_h)
grid::grid.newpage()
vp <- grid::viewport(
  x      = grid::unit(margin_left, "in"),
  y      = grid::unit(page_h - margin_top, "in"),
  width  = grid::unit(plot_w, "in"),
  height = grid::unit(plot_h, "in"),
  just   = c("left", "top"))
grid::pushViewport(vp)
grid::grid.draw(combined)
grid::popViewport()
dev.off()

message("Wrote: ", out_file)
