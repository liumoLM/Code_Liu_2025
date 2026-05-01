#!/usr/bin/env Rscript

# Build a "Rosetta stone" mapping table linking the indel classification
# columns Koh_476, Koh_89, COSMIC_83, and long_visual across the n_per_dir largest
# annotated indel VCFs from the FMH and PCAWG cohorts.

suppressPackageStartupMessages({
  library(argparser)
  library(dplyr)
  library(readr)
})

p <- arg_parser("Build rosetta stone mapping table from annotated indel VCFs")
p <- add_argument(p, "--n-per-dir", type = "integer", default = 300,
                  help = "Number of largest VCFs to use from each directory")
p <- add_argument(p, "--flank-5", type = "integer", default = 5,
                  help = "Number of characters to keep from the 5' flank in long_visual")
p <- add_argument(p, "--flank-3", type = "integer", default = 20,
                  help = "Number of characters to keep from the 3' flank in long_visual")
args <- parse_args(p)
n_per_dir  <- args$n_per_dir
flank_5    <- args$flank_5
flank_3    <- args$flank_3

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
  df$long_visual <- shorten_long_visual(df$long_visual, flank_5, flank_3)
  df[, keep_cols, drop = FALSE]
}

# long_visual is "<5'-flank> <indel-token> <3'-flank>". Keep the last flank_5
# chars of the 5' flank, the indel token, and the first flank_3 chars of the
# 3' flank.
shorten_long_visual <- function(x, flank_5 = 5, flank_3 = 20) {
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
    substr(pre, pmax(1, nchar(pre) - flank_5 + 1), nchar(pre)),
    mid,
    substr(post, 1, flank_3)
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
