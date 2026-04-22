# one_page_83_89_476.R
#
# Build a three-panel (ID83 / ID89 / ID476) stacked layout and draw it onto
# the current graphics device page inside a positioned viewport. The caller
# owns the device (e.g. cairo_pdf) and any grid.newpage() calls between
# pages.

one_page_83_89_476 <- function(
  sig_83,
  sig_89,
  sig_476,
  title_83 = NULL,
  title_89 = NULL,
  title_476 = NULL,
  num_peak_labels = 4,
  base_size = 9.5,
  page_h = 11,      # letter portrait height (in)
  plot_w = 6.5,     # plot-block width       (in)
  plot_h = 6.0,     # plot-block height      (in)
  margin_left = 1,  # from left edge         (in)
  margin_top = 1,   # from top edge          (in)
  gap_pt = 9        # ~1/8 in between panels (pt)
) {
  p83 <- mSigPlot::plot_ID83(
    sig_83,
    plot_title = title_83,
    num_peak_labels = num_peak_labels,
    base_size = base_size
  )
  p89 <- mSigPlot::plot_ID89(
    sig_89,
    plot_title = title_89,
    num_peak_labels = num_peak_labels,
    base_size = base_size
  )
  p476 <- mSigPlot::plot_ID476(
    sig_476,
    plot_title = title_476,
    num_peak_labels = num_peak_labels,
    base_size = base_size
  )

  p83 <- p83 + ggplot2::theme(
    plot.margin = ggplot2::margin(0, 0, gap_pt, 0, "pt")
  )
  p89 <- p89 + ggplot2::theme(
    plot.margin = ggplot2::margin(0, 0, gap_pt, 0, "pt")
  )
  p476 <- p476 + ggplot2::theme(
    plot.margin = ggplot2::margin(0, 0, 0, 0, "pt")
  )

  combined <- gridExtra::arrangeGrob(
    grobs = list(p83, p89, p476),
    ncol = 1,
    nrow = 3
  )

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
  invisible(NULL)
}
