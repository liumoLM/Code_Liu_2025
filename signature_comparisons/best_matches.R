library(gridGraphics)
library(gridExtra)

source("best_match_one_column.R")

#' Find best matches for signatures against reference signatures
#' Compares each column of our_sigs against combined reference signatures
#' Creates a PDF for each signature with the best matches
#' @param plotit A function that takes (vector, title) and returns a grob
#' @param our_sigs_path Path to the file containing signatures to match
#' @param ... One or more paths to reference signature files (will be cbind'd)
best_matches <- function(plotit, out_dir, our_sigs_path, ...) {
  ref_paths <- list(...)

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  if (length(ref_paths) < 1) {
    stop("At least one reference signature file path is required")
  }

  # Read our signatures
  our_sigs <- read.table(
    our_sigs_path,
    header = TRUE,
    sep = "\t",
    row.names = 1
  )

  # Read and combine reference signatures
  other_sigs_list <- lapply(ref_paths, function(path) {
    read.table(
      path,
      header = TRUE,
      sep = "\t",
      row.names = 1
    )
  })

  other_sigs <- do.call(cbind, other_sigs_list)

  # Convert to matrices
  other_sigs_mat <- as.matrix(other_sigs)
  our_sigs_mat <- as.matrix(our_sigs)
  results_list = list()

  # Process each column of our_sigs
  for (i in seq_len(ncol(our_sigs_mat))) {
    col_name <- colnames(our_sigs_mat)[i]
    A <- our_sigs_mat[, i, drop = FALSE]

    # Create PDF
    pdf_name <- sprintf("%s/%s_best_match.pdf", out_dir, col_name)
    cairo_pdf(pdf_name, width = 14, height = 12)

    result <- best_match_one_column(A, other_sigs_mat, plotit)

    dev.off()

    message(sprintf("Created %s", pdf_name))

    # Print results summary
    message(sprintf("Best matches for %s:", col_name))
    for (measure in names(result$results)) {
      r <- result$results[[measure]]
      message(sprintf("  %s: %s (%.4f)", measure, r$column_name, r$value))
    }

    # Build row for this signature
    row_data <- list(signature = col_name)
    for (measure in names(result$results)) {
      r <- result$results[[measure]]
      row_data[['ID']] <- r$column_name
      row_data[[measure]] <- r$value
    }
    results_list[[i]] <- tibble::as_tibble(row_data)
  }

  invisible(dplyr::bind_rows(results_list))
}
