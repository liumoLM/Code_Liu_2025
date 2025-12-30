library(ICAMS)
library(gridExtra)
library(mSigPlot)

source("code/wrap_ICAMS_plot_catalog.R")

#' Plot selected columns from a signature file to PDF
#' @param input_file_path Path to tab-separated input file
#' @param output_file_path Path for output PDF
#' @param identifiers Character vector of patterns to grep for in column names
#' @return Invisible NULL, creates PDF as side effect
plot_selected_columns <- function(
  input_file_path,
  output_file_path,
  identifiers
) {
  # Read input file
  data <- read.table(
    input_file_path,
    header = TRUE,
    sep = "\t",
    row.names = 1
  )

  n_rows <- nrow(data)

  # Find matching columns
  matched_cols <- unlist(lapply(identifiers, function(id) {
    grep(id, colnames(data), value = TRUE)
  }))
  matched_cols <- unique(matched_cols)

  if (length(matched_cols) == 0) {
    stop("No columns matched the provided identifiers")
  }

  # Determine plotting function based on row count
  if (n_rows == 83) {
    plot_fn <- function(vec, title) {
      catalog <- matrix(vec, ncol = 1)
      rownames(catalog) <- rownames(data)
      colnames(catalog) <- title
      catalog <- ICAMS::as.catalog(catalog, catalog.type = "counts.signature")
      wrap_ICAMS_plot_catalog(catalog, title)
    }
  } else if (n_rows == 89) {
    plot_fn <- function(vec, title) {
      plot_89(vec, plot_title = title)
    }
  } else if (n_rows == 476) {
    plot_fn <- function(vec, title) {
      plot_476(vec, plot_title = title)
    }
  } else {
    stop(sprintf(
      "Unsupported number of rows: %d. Expected 83, 89, or 476.",
      n_rows
    ))
  }

  # Create plots
  plot_list <- lapply(matched_cols, function(col_name) {
    vec <- as.numeric(data[, col_name])
    plot_fn(vec, col_name)
  })

  # Save to PDF - 1 column, 6 rows per page
  plots_per_page <- 6
  total_pages <- ceiling(length(plot_list) / plots_per_page)

  cairo_pdf(output_file_path, width = 10, height = 14, onefile = TRUE)

  for (page in seq_len(total_pages)) {
    start_idx <- (page - 1) * plots_per_page + 1
    end_idx <- min(page * plots_per_page, length(plot_list))
    plots_on_page <- plot_list[start_idx:end_idx]

    gridExtra::grid.arrange(grobs = plots_on_page, ncol = 1, nrow = 6)
  }

  dev.off()

  message(sprintf(
    "Created %s with %d plots",
    output_file_path,
    length(plot_list)
  ))
  invisible(NULL)
}
