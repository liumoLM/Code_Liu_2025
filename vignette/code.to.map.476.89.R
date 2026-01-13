ID476_ID89_mapping <- data.table::fread("./ID476_ID89_mapping.txt")

Convert_Indel476_to_Indel89 <- function(indel476.catalog){
  #row.names(indel476.catalog) <- indel476_catalogs$MutationType
  
  indel476.catalog$mut89_class <- ID476_ID89_mapping$indel89.class[match(
    row.names(indel476.catalog),ID476_ID89_mapping$indel476.class
  )]
  indel476.catalog$mut89_class[is.na(indel476.catalog$mut89_class)] <- "Complex"
  df_summary <- indel476.catalog %>%
    group_by(mut89_class) %>%
    summarise(across(where(is.numeric), ~sum(.x, na.rm = TRUE), .names = "{.col}_sum"))
  return(df_summary)
}