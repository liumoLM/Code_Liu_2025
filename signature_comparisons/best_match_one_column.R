library(gridExtra)
library(ggplot2)
library(indelsig.tools.lib)
library(philentropy)

#' Find best matching column in M for vector A using multiple similarity measures
#' @param A A column vector (or single-column matrix/data.frame)
#' @param M A matrix with columns to compare against
#' @param plotit A function that takes (vector, title) and returns a ggplot
#' @return A list with best matches and combined plot
best_match_one_column <- function(A, M, plotit) {
  # Ensure A is a numeric vector
  A_vec <- as.numeric(A)
  A_name <- if (!is.null(colnames(A))) {
    colnames(A)[1]
  } else {
    deparse(substitute(A))
  }

  # Normalize for probability-based measures
  A_norm <- A_vec / sum(A_vec)

  # Define similarity measures and whether higher is better
  measures <- list(
    cosine = list(method = "cosine", higher_better = TRUE),
    euclidean = list(method = "euclidean", higher_better = FALSE)
    # hellinger = list(method = "hellinger", higher_better = FALSE),
    # jaccard = list(method = "jaccard", higher_better = FALSE)
  )

  results <- list()
  plots <- list()

  # Plot A first
  plots[[1]] <- plotit(A_vec, A_name)

  for (measure_name in names(measures)) {
    measure <- measures[[measure_name]]
    best_idx <- NA
    best_val <- if (measure$higher_better) -Inf else Inf

    for (j in seq_len(ncol(M))) {
      M_col <- as.numeric(M[, j])

      # Normalize columns for comparison

      A_use <- A_norm
      M_use <- M_col / sum(M_col)

      # philentropy::distance expects a matrix with rows as distributions
      dist_mat <- rbind(A_use, M_use)

      val <- tryCatch(
        {
          suppressMessages(philentropy::distance(
            dist_mat,
            method = measure$method,
            test.na = FALSE
          )[1])
        },
        error = function(e) NA
      )

      if (!is.na(val)) {
        if (
          (measure$higher_better && val > best_val) ||
            (!measure$higher_better && val < best_val)
        ) {
          best_val <- val
          best_idx <- j
        }
      }
    }

    results[[measure_name]] <- list(
      index = best_idx,
      column_name = colnames(M)[best_idx],
      value = best_val
    )

    # Create plot for this match
    M_best <- as.numeric(M[, best_idx])
    title <- sprintf(
      "%s (%s: %.4f)",
      colnames(M)[best_idx],
      measure_name,
      best_val
    )
    plots[[length(plots) + 1]] <- plotit(M_best, title)
  }

  # Combine plots
  combined_plot <- gridExtra::grid.arrange(grobs = plots, ncol = 1)

  list(
    results = results,
    plots = plots,
    combined_plot = combined_plot
  )
}
