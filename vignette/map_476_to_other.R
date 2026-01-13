library(dplyr)

t476_to_89 <- function(t476) {
  mut_type_mapping <- data.table::fread("./ID476_ID89_mapping.txt")
  tmp <- read.delim(
    "../Manuscript_data/Liu_et_al_final_89_type_signatures.tsv",
    sep = '\t'
  )
  correct_row_order = tmp[, 1]
  rm(tmp)

  t476$mut89_class <- mut_type_mapping$indel89.class[match(
    row.names(t476),
    mut_type_mapping$indel476.class
  )]

  stopifnot(!is.na(t476$mut89_class))

  new89 <- t476 %>%
    group_by(mut89_class) %>%
    summarise(across(
      where(is.numeric),
      ~ sum(.x, na.rm = TRUE),
      .names = "{.col}_converted"
    )) %>%
    as.data.frame()

  stopifnot(length(symdiff(correct_row_order, new89[, 1])) == 0)

  row.names(new89) <- new89$mut89_class
  new89 <- new89[, -1, drop = FALSE]
  new89 = new89[correct_row_order, ]

  new89
}
