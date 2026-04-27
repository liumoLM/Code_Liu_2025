#!/usr/bin/env Rscript
# Reverse set difference: Koh_476 values present in
# rosetta_stone_476_89_cap9.csv but absent from indel476.class in
# ID476_ID89_mapping.txt.

suppressPackageStartupMessages({
  library(readr)
})

mapping <- read.delim(
  here::here("Manuscript_data", "ID476_ID89_mapping.txt"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

rosetta <- read.csv(
  here::here("some_sup_tables", "rosetta_stone_476_89_cap9.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

extra <- setdiff(rosetta$Koh_476, mapping$indel476.class)

out <- rosetta[match(extra, rosetta$Koh_476), , drop = FALSE]

write_tsv(out, here::here("some_sup_tables", "extra.tsv"))

message("Wrote ", nrow(out), " rows to some_sup_tables/extra.tsv")
