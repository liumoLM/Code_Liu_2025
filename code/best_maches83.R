library(ICAMS)
library(gridGraphics)
library(gridExtra)

source("code/wrap_ICAMS_plot_catalog.R")
source("code/best_matches.R")


#' Run best_matches for ID83 signatures against COSMIC and Jin 2024 references
#' @param out_dir Output directory for PDF files
#' @param ... At least 2 file paths: first is our_sigs, rest are reference sigs
best_matches83 <- function(out_dir, ...) {
  paths <- list(...)
  if (length(paths) < 2) {
    stop("At least 2 file paths required: our_sigs and at least one reference")
  }

  # Read row names from first file (our_sigs)
  our_sigs <- read.table(
    "data/type83_our_sigs.tsv",
    header = TRUE,
    sep = "\t",
    row.names = 1
  )
  rnames <- rownames(our_sigs)
  rm(our_sigs)

  # Define plotit function using ICAMS wrapper for ID83
  plotit <- function(vec, title) {
    catalog <- matrix(vec, ncol = 1)
    rownames(catalog) <- rnames
    colnames(catalog) <- title
    catalog <- ICAMS::as.catalog(catalog, catalog.type = "counts.signature")
    wrap_ICAMS_plot_catalog(catalog, title)
  }

  best_matches(
    plotit,
    out_dir,
    ...
  )
}

us_v_cosmic_83 = best_matches83(
  "83_us_v_cosmic",
  "data/type83_our_sigs.tsv",
  "data/COSMIC_v3.5_ID_GRCh37.txt"
)

us_v_jin_83 = best_matches83(
  "83_us_v_jin",
  "data/type83_our_sigs.tsv",
  "data/jin_2024_indel_sigs_sup_tab_1.tsv"
)

best_matches83(
  "83_jin_v_us",
  "data/jin_2024_indel_sigs_sup_tab_1.tsv",
  "data/type83_our_sigs.tsv"
)

best_matches83(
  "83_jin_v_cosmic",
  "data/jin_2024_indel_sigs_sup_tab_1.tsv",
  "data/COSMIC_v3.5_ID_GRCh37.txt"
)
# ...existing code...

us_vs_all <- dplyr::full_join(
  us_v_cosmic_83,
  us_v_jin_83,
  dplyr::join_by(signature)
) |>
  dplyr::mutate(
    max_cosine = pmax(cosine.x, cosine.y, na.rm = TRUE),
    max_cosine_id = dplyr::if_else(
      cosine.x >= cosine.y | is.na(cosine.y),
      ID.x,
      ID.y
    )
  ) |>
  dplyr::relocate(max_cosine_id, max_cosine, .after = signature) |>
  dplyr::arrange(dplyr::desc(max_cosine))
write.csv(us_vs_all, "signature_comparisons/us_vs_all_83.csv")
