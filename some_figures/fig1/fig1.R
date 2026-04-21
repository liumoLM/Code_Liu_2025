# fig1.R — 3-panel stacked indel-signature figure (ID83 / ID89 / ID476)
# Portrait letter PDF with editable text (cairo_pdf).
# Input: fig1_data.rds — a list with elements sig_83, sig_89, sig_476, each
# a single-column data frame of mutation-class counts/proportions.

library(ggplot2)

# ---- Params (edit these to tune the figure) -------------------------------

num_peaks <- 4

base_size <- 9.5
base_size_83 <- base_size
base_size_89 <- base_size
base_size_476 <- base_size

page_w <- 8.5 # letter portrait width  (in)
page_h <- 11 # letter portrait height (in)
plot_w <- 6.5 # plot-block width       (in)
plot_h <- 6.0 # plot-block height      (in) — ggrepel needs ~2 in per panel
margin_left <- 1 # from left edge         (in)
margin_top <- 1 # from top edge          (in)

# Locate the directory that contains this script so the .rds sibling
# resolves whether sourced interactively, via Rscript, or via source(chdir=).
this_script <- tryCatch(
  {
    ca <- commandArgs(trailingOnly = FALSE)
    m <- regmatches(ca, regexpr("(?<=--file=).+", ca, perl = TRUE))
    if (length(m)) {
      m[1]
    } else {
      # Walk every active frame looking for source()'s `ofile`. `sys.frame(1)`
      # alone breaks when source() is called from inside another function
      # (e.g. RStudio's Source button, or any wrapper).
      ofile <- NULL
      for (i in seq_len(sys.nframe())) {
        f <- sys.frame(i)
        if (!is.null(f$ofile) && is.character(f$ofile) && nzchar(f$ofile)) {
          ofile <- f$ofile
          break
        }
      }
      ofile
    }
  },
  error = function(e) NULL
)
script_dir <- if (!is.null(this_script) && nzchar(this_script)) {
  dirname(normalizePath(this_script))
} else {
  getwd()
}
data_file <- file.path(script_dir, "fig1_data.rds")
out_file <- file.path(script_dir, "fig1.pdf")

# ---- Read signatures ------------------------------------------------------

dat <- readRDS(data_file)

# ---- Build the three panels -----------------------------------------------

p83 <- mSigPlot::plot_ID83(
  dat$sig_83,
  num_peak_labels = num_peaks,
  base_size = base_size_83
)

p89 <- mSigPlot::plot_ID89(
  dat$sig_89,
  num_peak_labels = num_peaks,
  base_size = base_size_89
)

p476 <- mSigPlot::plot_ID476(
  dat$sig_476,
  num_peak_labels = num_peaks,
  base_size = base_size_476
)

# Strip outer plot margins so panels fill their 6.5" container horizontally.
# Add ~1/8 in (9 pt) of vertical space between panels via bottom margins
# on the upper two panels.
gap_pt <- 9 # 1/8 inch
p83 <- p83 +
  ggplot2::theme(plot.margin = ggplot2::margin(0, 0, gap_pt, 0, "pt"))
p89 <- p89 +
  ggplot2::theme(plot.margin = ggplot2::margin(0, 0, gap_pt, 0, "pt"))
p476 <- p476 + ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0, "pt"))

combined <- gridExtra::arrangeGrob(
  grobs = list(p83, p89, p476),
  ncol = 1,
  nrow = 3
)

# ---- Write PDF with plot block positioned in top printable area -----------

cairo_pdf(out_file, width = page_w, height = page_h)
grid::grid.newpage()
vp <- grid::viewport(
  x = grid::unit(margin_left, "in"),
  y = grid::unit(page_h - margin_top, "in"),
  width = grid::unit(plot_w, "in"),
  height = grid::unit(plot_h, "in"),
  just = c("left", "top")
)
grid::pushViewport(vp)
grid::grid.draw(combined)
grid::popViewport()
dev.off()

message("Wrote: ", out_file)
