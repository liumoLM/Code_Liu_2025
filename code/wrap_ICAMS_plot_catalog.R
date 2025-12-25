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

  ICAMS::PlotCatalog(catalog, xlabels = FALSE, upper = FALSE, ...)

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
