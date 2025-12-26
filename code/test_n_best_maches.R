source("code/best_n_matches.R")
source("code/wrap_ICAMS_plot_catalog.R")

our_sigs <- read.table(
  "data/type83_our_sigs.tsv",
  header = TRUE,
  sep = "\t",
  row.names = 1
)


our_sigs_mat <- as.matrix(our_sigs)

# Define plotit function using ICAMS wrapper
plotit_83 <- function(vec, title) {
  # Create a single-column catalog
  catalog <- matrix(vec, ncol = 1)
  rownames(catalog) <- rownames(our_sigs)
  colnames(catalog) <- title

  # Convert to ICAMS catalog
  catalog <- ICAMS::as.catalog(catalog, catalog.type = "counts.signature")

  wrap_ICAMS_plot_catalog(catalog, title)
}


# Get top 4 matches for each signature
results <- best_n_matches(
  plotit_83,
  "83_best_4_matches",
  "data/type83_our_sigs.tsv",
  4,
  "data/type83_spectra.tsv"
)
