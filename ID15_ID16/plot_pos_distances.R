#' Plot histogram of distances between successive genomic positions
#'
#' @param df Data frame with a POS column containing genomic positions
#' @param max_dist Maximum distance to include (default 10000)
#' @param bins Number of histogram bins (default 50)
#'
#' @return A ggplot2 object
plot_pos_distances <- function(df, max_dist = 10000, bins = 50) {
  library(ggplot2)

  # Calculate differences between successive POS values
  diffs <- diff(df$POS)

  # Filter to positive distances <= max_dist
  diffs_filtered <- diffs[diffs > 0 & diffs <= max_dist]

  ggplot(data.frame(dist = diffs_filtered), aes(x = dist)) +
    geom_histogram(bins = bins, fill = "steelblue", color = "white") +
    labs(
      title = paste0("Distribution of distances (<= ", max_dist, ")"),
      x = "Distance between successive positions",
      y = "Frequency"
    ) +
    theme_minimal()
}
