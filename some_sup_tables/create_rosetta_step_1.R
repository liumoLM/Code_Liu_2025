#!/usr/bin/env Rscript

n_per_dir <- 300

# Build a "Rosetta stone" mapping table linking the indel classification
# columns Koh_476, Koh_89, COSMIC_83, and long_visual across the n_per_dir largest
# annotated indel VCFs from the FMH and PCAWG cohorts.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

out_dir <- here::here("some_sup_tables")

vcf_dirs <- c(
  "/home/steve/MEGA/important_mut_sig_data/fmh-unfiltered_vcfs",
  "/home/steve/MEGA/important_mut_sig_data/pcawg_indel_vcfs"
)

keep_cols <- c(
  "Koh_476",
  "Koh_89",
  "COSMIC_83",
  "long_visual",
  "U_seq",
  "U_seq_count_in_indel_seq",
  "R",
  "mh",
  "unit",
  "spacer"
)

pick_largest <- function(dir, n) {
  files <- list.files(
    dir,
    pattern = "\\.annotated\\.indel\\.vcf\\.gz$",
    full.names = TRUE
  )
  if (length(files) == 0) {
    return(character(0))
  }
  sizes <- file.info(files)$size
  files[order(sizes, decreasing = TRUE)][seq_len(min(n, length(files)))]
}

read_one <- function(path) {
  message("Reading ", basename(path))
  df <- read.delim(
    gzfile(path),
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )
  missing <- setdiff(keep_cols, colnames(df))
  if (length(missing) > 0) {
    warning(path, " is missing columns: ", paste(missing, collapse = ", "))
    for (m in missing) {
      df[[m]] <- NA_character_
    }
  }
  df$long_visual <- shorten_long_visual(df$long_visual)
  df[, keep_cols, drop = FALSE]
}

# long_visual is "<5'-flank> <indel-token> <3'-flank>". Keep the last 5 chars
# of the 5' flank, the indel token, and the first 5 chars of the 3' flank.
shorten_long_visual <- function(x) {
  out <- x
  ok <- !is.na(x)
  parts <- strsplit(x[ok], " ", fixed = TRUE)
  ok_three <- vapply(parts, length, integer(1)) == 3L
  idx <- which(ok)[ok_three]
  good <- parts[ok_three]
  pre <- vapply(good, `[`, character(1), 1)
  mid <- vapply(good, `[`, character(1), 2)
  post <- vapply(good, `[`, character(1), 3)
  out[idx] <- paste(
    substr(pre, pmax(1, nchar(pre) - 4), nchar(pre)),
    mid,
    substr(post, 1, 5)
  )
  out
}

vcf_files <- unlist(lapply(vcf_dirs, pick_largest, n = n_per_dir))
message("Selected ", length(vcf_files), " VCF files")

big <- bind_rows(lapply(vcf_files, read_one))

write_csv(big, file.path(out_dir, "rosetta_stone_full.csv"))

unique_476_89 <- big |>
  dplyr::select(Koh_476, Koh_89) |>
  dplyr::distinct() |>
  dplyr::arrange(Koh_476, Koh_89)
write_csv(unique_476_89, file.path(out_dir, "rosetta_stone_476_89.csv"))

unique_476_89_83 <- big |>
  dplyr::select(Koh_476, Koh_89, COSMIC_83) |>
  dplyr::distinct() |>
  dplyr::arrange(Koh_476, Koh_89, COSMIC_83)
write_csv(unique_476_89_83, file.path(out_dir, "rosetta_stone_476_89_83.csv"))

message("Wrote:")
message(
  "  ",
  file.path(out_dir, "rosetta_stone_full.csv"),
  " (",
  nrow(big),
  " rows)"
)
message(
  "  ",
  file.path(out_dir, "rosetta_stone_476_89.csv"),
  " (",
  nrow(unique_476_89),
  " rows)"
)
message(
  "  ",
  file.path(out_dir, "rosetta_stone_476_89_83.csv"),
  " (",
  nrow(unique_476_89_83),
  " rows)"
)
