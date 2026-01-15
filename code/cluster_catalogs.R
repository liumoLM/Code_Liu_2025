library(ggplot2)
library(ggdendro)

#' Hierarchical clustering of signature catalogs using cosine distance
#'
#' @param output_file Path to the output PDF file
#' @param ... Named arguments where names are source labels and values are file paths
#'   to TSV files containing signature data (features as rows, samples as columns)
#'
#' @examples
#' # See code/run_cluster_catalogs.R for complete examples
cluster_catalogs <- function(output_file, ...) {
  file_args <- list(...)

  if (length(file_args) == 0) {
    stop("At least one file path must be provided")
  }

  # Check that all arguments are named

  arg_names <- names(file_args)
  if (is.null(arg_names) || any(arg_names == "")) {
    stop("All file arguments must be named (e.g., Source1 = 'path/to/file.tsv')")
  }

  # Read all files and track sources
  sig_list <- list()
  source_vec <- character(0)

  for (source_name in arg_names) {
    file_path <- file_args[[source_name]]
    sigs <- read.table(
      file_path,
      header = TRUE,
      sep = "\t",
      row.names = 1
    )
    sig_list[[source_name]] <- sigs
    source_vec <- c(source_vec, rep(source_name, ncol(sigs)))
  }

  # Combine all signatures and transpose (samples as rows, features as columns)
  combined <- do.call(cbind, sig_list)
  combined_t <- t(combined)
  names(source_vec) <- rownames(combined_t)

  # Normalize rows to unit vectors for cosine distance
  combined_norm <- combined_t / sqrt(rowSums(combined_t^2))

  # Compute cosine distance matrix (1 - cosine similarity)
  cosine_sim <- combined_norm %*% t(combined_norm)
  cosine_dist <- as.dist(1 - cosine_sim)

  # Hierarchical clustering (average linkage works well with cosine distance)
  hc <- hclust(cosine_dist, method = "average")

  # Extract dendrogram data
  dend_data <- dendro_data(hc)

  # Add source info to labels
  label_df <- dend_data$labels
  label_df$source <- source_vec[label_df$label]

  # Color-blind friendly palette (Okabe-Ito extended)
  cb_colors <- c(
    "#0072B2", "#D55E00", "#009E73", "#F0E442",
    "#CC79A7", "#56B4E9", "#E69F00", "#000000"
  )
  n_sources <- length(arg_names)
  if (n_sources > length(cb_colors)) {
    stop("Too many sources (max ", length(cb_colors), ")")
  }
  cb_palette <- setNames(cb_colors[seq_len(n_sources)], arg_names)

  # Plot
  p <- ggplot() +
    geom_segment(
      data = dend_data$segments,
      aes(x = x, y = y, xend = xend, yend = yend)
    ) +
    geom_text(
      data = label_df,
      aes(x = x, y = y - 0.01, label = label, color = source),
      hjust = 1,
      angle = 90,
      size = 3
    ) +
    geom_hline(yintercept = 0.1, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = cb_palette) +
    coord_cartesian(ylim = c(-0.15, NA), clip = "off") +
    labs(
      title = "Hierarchical Clustering of Signatures (Cosine Distance)",
      x = "",
      y = "Distance",
      color = "Source"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.margin = margin(t = 10, r = 10, b = 60, l = 10)
    )

  # Save to PDF (8.5 x 11 inch landscape)
  cairo_pdf(output_file, width = 11, height = 8.5)
  print(p)
  dev.off()

  message("Created ", output_file)

  invisible(list(
    hclust = hc,
    distance = cosine_dist,
    plot = p
  ))
}
