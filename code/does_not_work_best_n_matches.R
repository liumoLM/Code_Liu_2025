library(ICAMS)
library(gridGraphics)
library(gridExtra)
library(philentropy)

#' Find best N matches for signatures against reference signatures
#' Compares each column of our_sigs against combined reference signatures
#' Creates a PDF for each signature with the N best matches (cosine similarity)
#' @param plotit A function that takes (vector, title) and returns a plot/grob
#' @param out_dir Output directory for PDF files
#' @param our_sigs_path Path to the file containing signatures to match
#' @param num_best_matches Number of best matches to include in each plot
#' @param ... One or more paths to reference signature files (will be cbind'd)
best_n_matches <- function(
  plotit,
  out_dir,
  our_sigs_path,
  num_best_matches,
  ...
) {
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
  results_list <- list()

  # Process each column of our_sigs
  for (i in seq_len(ncol(our_sigs_mat))) {
    col_name <- colnames(our_sigs_mat)[i]
    A <- as.numeric(our_sigs_mat[, i])
    A_norm <- A / sum(A)

    # Compute cosine similarity to all reference columns
    cosine_vals <- numeric(ncol(other_sigs_mat))
    for (j in seq_len(ncol(other_sigs_mat))) {
      M_col <- as.numeric(other_sigs_mat[, j])
      M_norm <- M_col / sum(M_col)
      dist_mat <- rbind(A_norm, M_norm)
      cosine_vals[j] <- suppressMessages(
        philentropy::distance(dist_mat, method = "cosine", test.na = FALSE)[1]
      )
    }

    # Get indices of best matches (smallest cosine distance = highest similarity)
    n_matches <- min(num_best_matches, ncol(other_sigs_mat))
    best_indices <- order(cosine_vals)[1:n_matches]

    # Build plot list
    plot_list <- list()

    # First plot: the query signature
    plot_list[[1]] <- plotit(A, col_name)

    # Add best matches
    for (k in seq_along(best_indices)) {
      idx <- best_indices[k]
      match_name <- colnames(other_sigs_mat)[idx]
      match_vec <- as.numeric(other_sigs_mat[, idx])
      similarity <- 1 - cosine_vals[idx] # Convert distance to similarity
      title <- sprintf("%s (cosine: %.4f)", match_name, similarity)
      plot_list[[k + 1]] <- plotit(match_vec, title)
    }

    # Create PDF with 5 plots per page
    pdf_name <- sprintf(
      "%s/%s_best_%d_matches.pdf",
      out_dir,
      col_name,
      num_best_matches
    )
    cairo_pdf(pdf_name, width = 10, height = 14, onefile = TRUE)

    plots_per_page <- 5
    total_plots <- length(plot_list)
    total_pages <- ceiling(total_plots / plots_per_page)

    for (page in seq_len(total_pages)) {
      start_idx <- (page - 1) * plots_per_page + 1
      end_idx <- min(page * plots_per_page, total_plots)
      plots_on_page <- plot_list[start_idx:end_idx]

      gridExtra::grid.arrange(grobs = plots_on_page, ncol = 1, nrow = 5)
    }

    dev.off()

    message(sprintf("Created %s", pdf_name))

    # Build results row
    row_data <- list(signature = col_name)
    for (k in seq_along(best_indices)) {
      idx <- best_indices[k]
      row_data[[paste0("match_", k)]] <- colnames(other_sigs_mat)[idx]
      row_data[[paste0("cosine_", k)]] <- 1 - cosine_vals[idx]
    }
    results_list[[i]] <- tibble::as_tibble(row_data)
  }

  invisible(dplyr::bind_rows(results_list))
}
