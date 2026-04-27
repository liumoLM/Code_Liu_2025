#!/usr/bin/env Rscript
# Set difference: indel476.class values in ID476_ID89_mapping.txt that are
# not present in the Koh_476 column of rosetta_stone_476_89.csv.

suppressPackageStartupMessages({
  library(readr)
})

mapping <- read.delim(
  here::here("Manuscript_data", "ID476_ID89_mapping.txt"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

rosetta <- read.csv(
  here::here("some_sup_tables", "rosetta_stone_476_89.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

needed <- setdiff(mapping$indel476.class, rosetta$Koh_476)

out <- mapping[match(needed, mapping$indel476.class), , drop = FALSE]

write_tsv(out, here::here("some_sup_tables", "needed.tsv"))

message("Wrote ", nrow(out), " rows to some_sup_tables/needed.tsv")
