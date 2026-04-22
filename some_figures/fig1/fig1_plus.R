# fig1.R — 3-panel stacked indel-signature figure (ID83 / ID89 / ID476)
# Portrait letter PDF with editable text (cairo_pdf).
# Input: fig1_data.rds — a list with elements sig_83, sig_89, sig_476, each
# a single-column data frame of mutation-class counts/proportions.

library(ggplot2)

# ---- Params (edit these to tune the figure) -------------------------------

# When TRUE, draw peak-label arrows and write fig1_with_arrows.pdf.
# When FALSE, suppress all peak labels and write fig1.pdf.
plot_arrows <- FALSE

num_peaks <- if (plot_arrows) 4 else 0

base_size <- 9.5

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
out_file <- file.path(
  script_dir,
  "fig1.pdf"
)

# ---- Read signatures ------------------------------------------------------

dat <- readRDS(data_file)

# ---- Write PDF with plot block positioned in top printable area -----------

source(file.path(script_dir, "one_page_83_89_476.R"))

mytxt <- function(label, x, y) {
  grid::grid.text(
    label = label,
    x = unit(x, "pt"),
    y = unit(y, "pt"),
    gp = gpar(fontsize = base_size, lineheight = 1)
  )
}

myseg <- function(x0, y0, x1, y1) {
  grid::grid.segments(
    x0 = unit(x0, "pt"),
    y0 = unit(y0, "pt"),
    x1 = unit(x1, "pt"),
    y1 = unit(y1, "pt")
  )
}

cairo_pdf(out_file, width = page_w, height = page_h)
one_page_83_89_476(
  title_83 = "83-Type Classification",
  title_89 = "89-Type Classification",
  title_476 = "476-Type Classification",
  sig_83 = dat$sig_83,
  sig_89 = dat$sig_89,
  sig_476 = dat$sig_476,
  num_peak_labels = num_peaks,
  base_size = base_size,
  page_h = page_h,
  plot_w = plot_w,
  plot_h = plot_h,
  margin_left = margin_left,
  margin_top = margin_top
)

# ID 83 =============

mytxt("VTV\u2192VV", 131, 681)
mytxt("VTTV\u2192VTV", 170, 665)
myseg(141, 635, 162, 659)
mytxt("e.g. TATA\u2192TA\ne.g. TGTG\u2192TG", 240, 640)
mytxt("2 base deletion\nwith microhomology\ne.g. GAG\u2192G", 484, 642)

# ID 89 =============

yA <- 550
myseg(158, yA, 168, yA)
mytxt(
  "AT{1,4}A\u2192AT{0,3}A, i.e. ATA\u2192AA, ATTA\u2192ATA, ATTTA\u2192ATTA, or ATTTTA\u2192ATTTA",
  345,
  yA
)

myseg(169.4, 529.5, 178.5, 535)
mytxt("CT{1,4}A\u2192CT{0,3}A", 225, 535)

yGT <- 520
myseg(189, yGT, 194, yGT)
mytxt("GT{1,4}A\u2192GT{0,3}A", 240, yGT)

yM <- 526
myseg(443.3, yM, 453, yM)
mytxt("Many possibilities", 405, yM)
mytxt("Many\npossibilities", 507, 530)

# ID 476 =============

mytxt("ATA\u2192AA", 186, 404)
mytxt("CTA\u2192CA", 218, 385)
mytxt("GTA\u2192GA", 260, 359)
mytxt("e.g. TATA\u2192TA\ne.g. TGTG\u2192TG", 416, 398)
# mytxt("e.g. TATATA\n\u2192TATA", 460, 345)
mytxt("2 base deletion\nwith microhomology\ne.g. GAG\u2192G", 503, 381)

dev.off()

message("Wrote: ", out_file)
