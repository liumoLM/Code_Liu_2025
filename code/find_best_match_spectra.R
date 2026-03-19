library(lsa)

#' Find the best-matching spectrum for each signature from a catalog
#'
#' @param signatures Matrix with mutation types as rows and signatures as columns
#' @param catalog_path Path to a TSV catalog file (mutation types as rows, samples as columns)
#' @return Subset of the catalog containing only the unique best-match columns (raw counts)
find_best_match_spectra <- function(signatures, catalog_path) {
  catalog <- read.table(
    catalog_path,
    header = TRUE,
    sep = "\t",
    row.names = 1,
    check.names = FALSE
  )
  cat_mat <- as.matrix(catalog)
  best_samples <- character(0)
  for (i in seq_len(ncol(signatures))) {
    sig_vec <- signatures[, i]
    cos_sims <- apply(cat_mat, 2, function(x) lsa::cosine(sig_vec, x))
    best_samples <- c(best_samples, names(which.max(cos_sims)))
  }
  unique_samples <- unique(best_samples)
  catalog[, unique_samples, drop = FALSE]
}

#' Find the names of the top N best-matching samples for each signature
#'
#' @param signatures Matrix with mutation types as rows and signatures as columns
#' @param catalog_path Path to a TSV catalog file (mutation types as rows, samples as columns)
#' @param n Number of top matches to return per signature
#' @return Character matrix with one row per signature and n columns
#'   (named BestMatch_1, BestMatch_2, ...). Row names are the signature column names.
find_top_n_match_names <- function(signatures, catalog_path, n = 3) {
  catalog <- read.table(
    catalog_path,
    header = TRUE,
    sep = "\t",
    row.names = 1,
    check.names = FALSE
  )
  cat_mat <- as.matrix(catalog)
  result <- matrix(
    NA_character_,
    nrow = ncol(signatures),
    ncol = n,
    dimnames = list(colnames(signatures), paste0("BestMatch_", seq_len(n)))
  )
  for (i in seq_len(ncol(signatures))) {
    sig_vec <- signatures[, i]
    cos_sims <- apply(cat_mat, 2, function(x) lsa::cosine(sig_vec, x)[1, 1])
    top_n <- names(sort(cos_sims, decreasing = TRUE)[seq_len(min(n, length(cos_sims)))])
    result[i, seq_along(top_n)] <- top_n
  }
  result
}
