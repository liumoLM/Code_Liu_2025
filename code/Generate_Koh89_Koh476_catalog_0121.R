##Input

## this function directly takes output from ICAMS::AnnovateIDVcf.
#' @param muts_list ICAMS annotated ID vcf
#' @param sample_col  the column containing the sample name
GenerateKoh89CatalogfromAnnotateVcf <- function(muts_list, sample_col) {
  template89 <- data.frame(IndelType = ICAMS::catalog.row.order$ID89)
  muts_list <- as.data.frame(muts_list)
  indel_catalogue <- data.frame(table(
    muts_list[, sample_col],
    muts_list$Koh_89
  ))
  names(indel_catalogue) <- c("Sample", "IndelType", "freq")
  indel_catalogue <- reshape2::dcast(
    indel_catalogue,
    IndelType ~
      Sample,
    value.var = "freq"
  )
  indel_catalogue <- merge(
    template89,
    indel_catalogue,
    by = "IndelType",
    all.x = T
  )
  indel_catalogue[is.na(indel_catalogue)] <- 0
  rownames(indel_catalogue) <- indel_catalogue[, "IndelType"]
  return(indel_catalogue[template89$IndelType, -1, drop = FALSE])
}

## this function directly takes output from ICAMS::AnnovateIDVcf.
#' @param muts_list ICAMS annotated ID vcf
#' @param sample_col  the column containing the sample name

GenerateKoh476CatalogfromAnnotateVcf <- function(muts_list, sample_col) {
  template476 <- data.frame(IndelType = ICAMS::catalog.row.order$ID476)
  muts_list <- as.data.frame(muts_list)
  indel_catalogue <- data.frame(table(
    muts_list[, sample_col],
    muts_list$Koh_476
  ))
  names(indel_catalogue) <- c("Sample", "IndelType", "freq")
  indel_catalogue <- reshape2::dcast(
    indel_catalogue,
    IndelType ~
      Sample,
    value.var = "freq"
  )
  indel_catalogue <- merge(
    template476,
    indel_catalogue,
    by = "IndelType",
    all.x = T
  )
  indel_catalogue[is.na(indel_catalogue)] <- 0
  rownames(indel_catalogue) <- indel_catalogue[, "IndelType"]
  return(indel_catalogue[template476$IndelType, -1, drop = FALSE])
}
