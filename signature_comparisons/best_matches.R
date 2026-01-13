library(gridGraphics)
library(gridExtra)

source("best_match_one_column.R")

guess_plotit = function(our_sigs_path) {
  our_sigs <- read.table(
    our_sigs_path,
    header = TRUE,
    sep = "\t",
    row.names = 1
  )
  if (nrow(our_sigs) == 89) {
    mSigPlot::plot_89
  } else if (nrow(our_sigs) == 83) {
    mSigPlot::plot_83
  } else if (nrow(our_sigs) == 476) {
    mSigPlot::plot_476
  } else {
    stop("unexpected number of row ", nrow(our_sigs))
  }
}

#' Find best matches for signatures against other signatures
#' or against spectra.
#' Compares each column of our_sigs against combined reference signatures
#' Creates a PDF for each signature with the best matches
#' @param plotit A function that takes (vector, title) and returns a grob
#' @param our_sigs_path Path to the file containing signatures to match
#' @param ... One or more paths to reference signature files (will be cbind'd)
best_matches <- function(
  plotit = guess_plotit(our_sigs_path),
  out_dir,
  min_num_mutations = 50,
  our_sigs_path,
  comp_path
) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  our_sigs <- read.table(
    our_sigs_path,
    header = TRUE,
    sep = "\t",
    row.names = 1
  )

  other_sigs <- read.table(
    comp_path,
    header = TRUE,
    sep = "\t",
    row.names = 1
  )

  # Convert to matrices
  other_sigs_mat <- as.matrix(other_sigs)
  our_sigs_mat <- as.matrix(our_sigs)
  results_list = list()

  # Process each column of our_sigs
  for (i in seq_len(ncol(our_sigs_mat))) {
    col_name <- colnames(our_sigs_mat)[i]
    A <- our_sigs_mat[, i, drop = FALSE]

    # Create PDF
    pdf_name <- sprintf(
      "%s/%s_%d_best_match.pdf",
      out_dir,
      col_name,
      nrow(our_sigs_mat)
    )
    cairo_pdf(pdf_name, width = 14, height = 12)

    result <- best_match_one_column(
      A,
      other_sigs_mat,
      min_num_mutations = min_num_mutations,
      plotit = plotit
    )

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

runit = FALSE
if (runit) {
  outdir0 = "best_spectra_matches"

  best_matches(
    out_dir = file.path(outdir0, "plots"),
    min_num_mutations = 50,
    our_sigs_path = "../Manuscript_data/Liu_et_al_final_89_type_signatures.tsv",
    comp_path = "../Manuscript_data/Liu_et_al_89_type_spectra.tsv"
  ) -> t89

  write.csv(t89, file.path(outdir0, "89_50_mutations.csv"))

  best_matches(
    out_dir = file.path(outdir0, "plots"),
    min_num_mutations = 50,
    our_sigs_path = "../Manuscript_data/Liu_et_al_final_476_type_signatures.tsv",
    comp_path = "../Manuscript_data/Liu_et_al_476_type_spectra.tsv"
  ) -> t476

  write.csv(t89, file.path(outdir0, "476_50_mutations.csv"))

  best_matches(
    out_dir = file.path(outdir0, "plots"),
    min_num_mutations = 50,
    our_sigs_path = "../Manuscript_data/Liu_et_al_final_83_type_signatures.tsv",
    comp_path = "../Manuscript_data/Liu_et_al_83_type_spectra.tsv"
  ) -> t83

  write.csv(t89, file.path(out_dir0, "83_50_mutations.csv"))
}
