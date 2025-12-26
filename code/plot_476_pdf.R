library(gridExtra)
library(Cairo)

plot_476_pdf <- function(
  Koh476_catalog,
  filename,
  num_labels = 4,
  simplify_labels = FALSE,
  label_threshold_denominator = 7
) {
  plot_list <- lapply(1:ncol(Koh476_catalog), function(i) {
    plot_476_v3(
      Koh476.catalog = Koh476_catalog[, i],
      text_size = 3,
      plot_title = colnames(Koh476_catalog)[i],
      num_labels = num_labels,
      simplify_labels = simplify_labels,
      label_threshold_denominator = label_threshold_denominator
    )
  })
  plots_per_page <- 5

  # Total number of pages
  total_pages <- ceiling(length(plot_list) / plots_per_page)

  # Open a PDF device
  cairo_pdf(filename, width = 8.2677, height = 14.61613)

  # Loop through pages and save each group of plots to the same PDF
  for (page in 1:total_pages) {
    # Get the plots for the current page
    start_index <- (page - 1) * plots_per_page + 1
    end_index <- min(page * plots_per_page, length(plot_list))
    plots_on_page <- plot_list[start_index:end_index]

    # Arrange the plots in a grid (2 rows x 4 columns)
    do.call(grid.arrange, c(plots_on_page, nrow = 5, ncol = 1))
  }

  # Close the PDF device
  dev.off()
}
