annot_vcf_to_476_catalog <- function(
  annot_vcf,
  sample_id = "no_sample_id_provided",
  FILTER_PASS = FALSE
) {
  domessage = TRUE
  if (colnames(annot_vcf)[1] == "CHROM") {
    colnames(annot_vcf)[1] <- "CHROM"
  }

  if (domessage) {
    message("initial annot_vcf rows = ", nrow(annot_vcf))
  }

  if (FILTER_PASS) {
    annot_vcf %>%
      dplyr::filter(FILTER == "PASS") -> annot_vcf
  }

  if (domessage) {
    message("num PASS rows = ", nrow(annot_vcf))
  }

  annot_vcf %>%
    dplyr::mutate(pos_id = paste0(CHROM, "-", POS)) -> vcf_with_pos_id

  if (domessage) {
    message("Checking for dup rows with different ALT values")
  }
  vcf_with_pos_id %>%
    dplyr::group_by(pos_id) %>%
    dplyr::filter(n_distinct(ALT) > 1) %>%
    select(pos_id, REF, ALT) -> multiple_alts

  if (nrow(multiple_alts) > 0) {
    warning(
      "Differences in 'ALT'; only 1 ALT value chosen arbitrarily at the following positions: ",
      paste(capture.output(print(multiple_alts)), collapse = '\n')
    )
  }

  vcf_with_pos_id %>% dplyr::distinct(pos_id, .keep_all = TRUE) -> cleaner_vcf
  if (domessage) {
    message("num PASS && unique rows = ", nrow(cleaner_vcf))
  }

  cleaner_vcf %>%
    dplyr::mutate(
      Koh_476 = if_else(
        R >= 9 &
          stringr::str_detect(
            Koh_476,
            "Del\\(T\\)|Del\\(C\\)|Ins\\(C\\)|Ins\\(T\\)"
          ),
        stringr::str_replace(Koh_476, "R\\d+", "R(9,)"),
        Koh_476
      )
    ) -> update_type_strings

  if (domessage) {
    message(
      "num PASS && unique rows after updating type string = ",
      nrow(update_type_strings)
    )
  }

  update_type_strings %>%
    dplyr::count(Koh_476) -> compacted_vcf

  if (domessage) {
    message("num PASS && unique mutations = ", sum(compacted_vcf$n))
  }

  data.table::data.table(Koh_476 = ICAMS::catalog.row.order$ID476) %>%
    dplyr::left_join(compacted_vcf) %>%
    mutate(n = if_else(is.na(n), 0L, n)) -> almost
  # browser()
  rownames(almost) <- pull(almost, Koh_476)
  colnames(almost)[2] <- sample_id

  almost[, -1, drop = FALSE]
}
