check_exemplars = function() {
  library(dplyr)

  table1 <- read.delim(
    "../Manuscript_data/table_1_2026_01_02.csv",
    sep = "\t",
    stringsAsFactors = FALSE
  ) %>%
    select(signature, Exemplar, maybe_update) %>%
    dplyr::filter(!signature %in% c("InsDel15", "InsDel16"))

  min_100 <- read.delim(
    "us_v_spectra_89_min_100.tsv",
    sep = "\t",
    stringsAsFactors = FALSE,
    row.names = 1
  ) %>%
    mutate(ID = gsub("..", "::", ID, fixed = TRUE))
  colnames(min_100) <- paste0(colnames(min_100), "_min_100")
  rownames(min_100) <- NULL

  no_min <- read.delim(
    "us_v_spectra_89.tsv",
    sep = "\t",
    stringsAsFactors = FALSE,
    row.names = 1
  ) %>%
    mutate(ID = gsub("..", "::", ID, fixed = TRUE))
  colnames(no_min) <- paste0(colnames(no_min), "_no_min")
  rownames(no_min) <- NULL

  result <- table1 %>%
    left_join(
      min_100,
      by = dplyr::join_by(signature == signature_min_100)
    ) %>%
    left_join(no_min, by = dplyr::join_by(signature == signature_no_min))

  return(filter(result, Exemplar != ID_min_100))
}
