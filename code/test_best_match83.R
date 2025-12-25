library(ICAMS)
library(gridGraphics)
library(gridExtra)

source("code/wrap_ICAMS_plot_catalog.R")
source("code/best_match_one_column.R")

#' Test best_match using ID83 signatures
#' Compares each column of type83_our_sigs against COSMIC ID signatures
#' Creates a PDF for each signature with the best matches
test_best_match83 <- function() {
  # Read the data files
  cosmic <- read.table(
    "data/COSMIC_v3.5_ID_GRCh37.txt",
    header = TRUE,
    sep = "\t",
    row.names = 1
  )

  our_sigs <- read.table(
    "data/type83_our_sigs.tsv",
    header = TRUE,
    sep = "\t",
    row.names = 1
  )

  # Convert to matrices
  cosmic_mat <- as.matrix(cosmic)
  our_sigs_mat <- as.matrix(our_sigs)

  # Define plotit function using ICAMS wrapper
  plotit <- function(vec, title) {
    # Create a single-column catalog
    catalog <- matrix(vec, ncol = 1)
    rownames(catalog) <- rownames(our_sigs)
    colnames(catalog) <- title

    # Convert to ICAMS catalog
    catalog <- ICAMS::as.catalog(catalog, catalog.type = "counts.signature")

    wrap_ICAMS_plot_catalog(catalog, title)
  }

  # Process each column of our_sigs
  for (i in 1:2) {
    col_name <- colnames(our_sigs_mat)[i]
    A <- our_sigs_mat[, i, drop = FALSE]

    # Create PDF
    pdf_name <- sprintf("testoutput/%s_best_match.pdf", col_name)
    pdf(pdf_name, width = 14, height = 12)

    result <- best_match_one_column(A, cosmic_mat, plotit)

    dev.off()

    message(sprintf("Created %s", pdf_name))

    # Print results summary
    message(sprintf("Best matches for %s:", col_name))
    for (measure in names(result$results)) {
      r <- result$results[[measure]]
      message(sprintf("  %s: %s (%.4f)", measure, r$column_name, r$value))
    }
  }
}

test_best_match83()
