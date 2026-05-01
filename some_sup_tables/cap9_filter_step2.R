#!/usr/bin/env Rscript
# Drop Koh_476 rows whose repeat count is >= 10 (i.e. the "uncapped" rows
# that have no counterpart in the cap-9 ID476_ID89 mapping). Write the
# filtered tables alongside the originals with "_cap9" before the suffix.

suppressPackageStartupMessages({
  library(data.table)
})

# Returns TRUE for Koh_476 strings like "X[Del(C):Rnn]Y", "X[Del(T):Rnn]Y",
# "X[Ins(C):Rnn]Y", or "X[Ins(T):Rnn]Y" where nn is a repeat count of 10+.
is_uncapped_del_repeat <- function(x) {
  grepl("^[ACGT]\\[(Del|Ins)\\([CT]\\):R1[0-9]+\\][ACGT]$", x)
}

in_dir <- here::here("some_sup_tables")
files <- c(
  "rosetta_stone_full.csv",
  "rosetta_stone_476_89.csv",
  "rosetta_stone_476_89_83.csv"
)

for (f in files) {
  in_path <- file.path(in_dir, f)
  out_path <- file.path(in_dir, sub("\\.csv$", "_cap9.csv", f))
  message("Reading ", f)
  dt <- fread(in_path, showProgress = FALSE)
  drop <- is_uncapped_del_repeat(dt$Koh_476)
  message("  dropping ", sum(drop), " of ", nrow(dt), " rows")
  dt <- dt[!drop]
  # Rewrite "(Ins|Del)(C|T):R9]" as "(Ins|Del)(C|T):R(9,)]" so the Koh_476
  # values match the ID476_ID89 mapping convention for the open-ended bin.
  dt$Koh_476 <- sub(
    "((?:Ins|Del)\\([CT]\\)):R9\\]",
    "\\1:R(9,)]",
    dt$Koh_476,
    perl = TRUE
  )
  if (f != "rosetta_stone_full.csv") {
    dt <- unique(dt)
  }
  fwrite(dt, out_path)
  message("  wrote ", out_path, " (", nrow(dt), " rows)")
}
