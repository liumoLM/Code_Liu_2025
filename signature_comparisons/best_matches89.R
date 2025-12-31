library(gridExtra)
library(mSigPlot)

source("code/best_matches.R")

#' Run best_matches for ID89 signatures
#' @param out_dir Output directory for PDF files
#' @param ... At least 2 file paths: first is our_sigs, rest are reference sigs
best_matches89 <- function(out_dir, ...) {
  paths <- list(...)
  if (length(paths) < 2) {
    stop("At least 2 file paths required: our_sigs and at least one reference")
  }

  # Define plotit function using PlotKoh89Catalog directly
  plotit <- function(vec, title) {
    plot_89(vec, plot_title = title)
  }

  best_matches(
    plotit,
    out_dir,
    ...
  )
}


us_v_koh_89 = best_matches89(
  "89_us_v_koh",
  "data/type89_our_sigs.tsv",
  "data/type89_koh_sigs.tsv"
) |>
  dplyr::arrange(desc(cosine))
write.csv(
  us_v_koh_89,
  "signature_comparisons/us_v_koh_89.csv",
  row.names = FALSE
)

koh_v_us_89 =
  best_matches89(
    "89_koh_v_us",
    "data/type89_koh_sigs.tsv",
    "data/type89_our_sigs.tsv"
  ) |>
  dplyr::arrange(desc(cosine))
write.csv(
  koh_v_us_89,
  "signature_comparisons/koh_v_us_89.csv",
  row.names = FALSE
)

us_v_spectra89 = best_matches89(
  "data/us_v_spectra89",
  "data/type89_liu_et_al_sigs.tsv",
  "data/type89_spectra.tsv"
)
write.table(us_v_spectra89, file = "data/us_v_spectra_80.tsv", sep = '\t')


koh_v_spectra89 = best_matches89(
  out_dir = "koh_v_spectra89",
  "../Manuscript_data/Koh_signatures.tsv",
  "../Manuscript_data/Liu_et_al_89_type_spectra.tsv"
)
write.table(koh_v_spectra89, file = "koh_v_spectra_80.tsv", sep = '\t')
