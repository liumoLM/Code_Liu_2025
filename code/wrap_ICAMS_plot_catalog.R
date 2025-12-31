library(ICAMS)
library(gridGraphics)

#' Wrap ICAMS::PlotCatalog to return a grob for use with grid.arrange
#' @param catalog A catalog object for ICAMS::PlotCatalog
#' @param ... Additional arguments passed to ICAMS::PlotCatalog
#' @return A grob that can be combined with ggplot objects via grid.arrange
wrap_ICAMS_plot_catalog <- function(catalog, title, ...) {
  stopifnot(ncol(catalog) == 1)
  colnames(catalog) <- title

  # Open a null PDF device
  pdf(NULL)
  dev.control(displaylist = "enable")

  # Draw the base R plot
  catalog.type = ifelse(sum(catalog) < 1.1, "counts.signature", "counts")

  c2 = ICAMS::as.catalog(
    catalog,
    catalog.type = catalog.type
  )

  ICAMS::PlotCatalog(c2, xlabels = TRUE, upper = FALSE, ...)

  # Record the base plot
  recorded_plot <- recordPlot()
  dev.off()
  # Replay and convert to grid
  pdf(NULL)
  dev.control(displaylist = "enable")
  replayPlot(recorded_plot)
  gridGraphics::grid.echo()
  grob <- grid::grid.grab()
  dev.off()

  grob
}

fix_cosmic_id = function(x) {
  stringi::stri_paste(
    stri_sub(x, 3, 8),
    stri_sub(x, 1, 2),
    stri_sub(x, 9, -1)
  ) |>
    stringi::stri_replace_all_fixed(":5", ":5+") |>
    stringi::stri_replace_all_fixed(":R:", ":repeats:") |>
    stringi::stri_replace_all_fixed(":M:", ":MH:") |>
    stringi::stri_replace_all_fixed("Del", "DEL") |>
    stringi::stri_replace_all_fixed("Ins", "INS")
}
