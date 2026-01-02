library(gridExtra)
library(ggplot2)
library(philentropy)

#' Find exemplar spectra for a given signature
#' @param sig_path Path to .tsv file with signatures
#' @param sig_col Column name of the signature to compare
#' @param spectra_path Path to .tsv file with mutational spectra
#' @param cosine_cutoff Minimum cosine similarity to include in results
#' @param out_pdf Path to output PDF file (default: sig_col_exemplars.pdf)
#' @param num_exemplars Number of top exemplars to plot (default: 10)
#' @param min_mutations Minimum total mutations to include spectrum (default: 0)
#' @return A data frame with cosine similarities for all spectra above cutoff
#' @export
exemplars_one_sig <- function(
  sig_path,
  sig_col,
  spectra_path,
  cosine_cutoff = 0.8,
  out_pdf = NULL,
  num_exemplars = 10,
  min_mutations = 0
) {
  # Read signature and spectra files
  sigs <- read.table(sig_path, sep = "\t", header = TRUE, row.names = 1)
  spectra <- read.table(spectra_path, sep = "\t", header = TRUE, row.names = 1)

  # Check that sig_col exists
  if (!sig_col %in% colnames(sigs)) {
    stop(paste("Column", sig_col, "not found in signature file"))
  }

  # Check row counts match
  if (nrow(sigs) != nrow(spectra)) {
    stop(paste(
      "Row counts don't match: signatures have",
      nrow(sigs),
      "rows, spectra have",
      nrow(spectra),
      "rows"
    ))
  }

  # Select plotting function based on number of rows
  nrows <- nrow(sigs)
  if (nrows == 83) {
    plotit <- function(vec, title) plot_83(vec, plot_title = title)
  } else if (nrows == 89) {
    plotit <- function(vec, title) plot_89(vec, plot_title = title)
  } else if (nrows == 476) {
    plotit <- function(vec, title) plot_476(vec, plot_title = title)
  } else {
    stop(paste(
      "Unsupported number of rows:",
      nrows,
      "(expected 83, 89, or 476)"
    ))
  }

  # Get signature vector and normalize
  sig_vec <- as.numeric(sigs[, sig_col])
  sig_norm <- sig_vec / sum(sig_vec)

  # Compute cosine similarity to each spectrum
  cosine_results <- data.frame(
    spectrum = character(),
    cosine = numeric(),
    total_mutations = numeric(),
    stringsAsFactors = FALSE
  )

  for (j in seq_len(ncol(spectra))) {
    spec_vec <- as.numeric(spectra[, j])
    total_muts <- sum(spec_vec)

    # Skip spectra with fewer than min_mutations
    if (total_muts < min_mutations) {
      next
    }

    spec_norm <- spec_vec / total_muts

    dist_mat <- rbind(sig_norm, spec_norm)

    cosine_val <- tryCatch(
      {
        suppressMessages(philentropy::distance(
          dist_mat,
          method = "cosine",
          test.na = FALSE
        )[1])
      },
      error = function(e) NA
    )

    if (!is.na(cosine_val)) {
      cosine_results <- rbind(
        cosine_results,
        data.frame(
          spectrum = colnames(spectra)[j],
          cosine = cosine_val,
          total_mutations = total_muts,
          stringsAsFactors = FALSE
        )
      )
    }
  }

  # Sort by cosine similarity (descending)
  cosine_results <- cosine_results[order(-cosine_results$cosine), ]

  # Filter by cutoff
  above_cutoff <- cosine_results[cosine_results$cosine >= cosine_cutoff, ]

  # Select top exemplars for plotting
  top_exemplars <- head(cosine_results, num_exemplars)

  # Set output PDF path
  if (is.null(out_pdf)) {
    out_pdf <- paste0(sig_col, "_exemplars.pdf")
  }

  # Create plots
  plots <- list()

  # Plot signature first
  sig_df <- sigs[, sig_col, drop = FALSE]
  plots[[1]] <- plotit(sig_df, paste("Signature:", sig_col))

  # Plot top exemplars
  for (i in seq_len(nrow(top_exemplars))) {
    spec_name <- top_exemplars$spectrum[i]
    cosine_val <- top_exemplars$cosine[i]
    total_muts <- top_exemplars$total_mutations[i]

    spec_df <- spectra[, spec_name, drop = FALSE]
    title <- sprintf(
      "%s (cosine: %.4f, n=%d)",
      spec_name,
      cosine_val,
      total_muts
    )

    plots[[length(plots) + 1]] <- plotit(spec_df, title)
  }

  # Save to PDF
  pdf_height <- 3 * length(plots)
  pdf_width <- 12

  cairo_pdf(out_pdf, width = pdf_width, height = pdf_height)
  combined_plot <- gridExtra::grid.arrange(grobs = plots, ncol = 1)
  dev.off()

  message(paste("Saved plot to", out_pdf))
  message(paste(
    "Found",
    nrow(above_cutoff),
    "spectra above cosine cutoff of",
    cosine_cutoff
  ))

  # Return results
  list(
    all_results = cosine_results,
    above_cutoff = above_cutoff,
    top_exemplars = top_exemplars
  )
}

exemplars_one_sig(
  "../Manuscript_data/COSMIC_v3.5_ID_GRCh37_signatures.tsv",
  "ID10",
  "../Manuscript_data/Liu_et_al_83_type_spectra.tsv",
  cosine_cutoff = 0.5,
  out_pdf = "exemplars_for_ID10.pdf",
  min_mutations = 50
)

exemplars_one_sig(
  "../Manuscript_data/COSMIC_v3.5_ID_GRCh37_signatures.tsv",
  "ID10",
  "../Manuscript_data/Liu_et_al_83_type_spectra.tsv",
  cosine_cutoff = 0.5,
  out_pdf = "exemplars_for_ID10_min_20.pdf",
  min_mutations = 20
)


exemplars_one_sig(
  "../Manuscript_data/Liu_et_al_final_83_type_signatures.tsv",
  "C_ID10",
  "../Manuscript_data/Liu_et_al_83_type_spectra.tsv",
  cosine_cutoff = 0.5,
  out_pdf = "exemplars_for_C_ID10.pdf",
  min_mutations = 50
)


exemplars_one_sig(
  "../Manuscript_data/Liu_et_al_final_83_type_signatures.tsv",
  "C_ID7",
  "../Manuscript_data/Liu_et_al_83_type_spectra.tsv",
  cosine_cutoff = 0.5,
  out_pdf = "exemplars_for_C_ID7_min_50.pdf",
  min_mutations = 50
)

exemplars_one_sig(
  "../Manuscript_data/COSMIC_v3.5_ID_GRCh37_signatures.tsv",
  "ID7",
  "../Manuscript_data/Liu_et_al_83_type_spectra.tsv",
  cosine_cutoff = 0.5,
  out_pdf = "exemplars_for_ID7_min_50.pdf",
  min_mutations = 50
)
