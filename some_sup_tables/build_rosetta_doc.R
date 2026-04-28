#!/usr/bin/env Rscript
# Build a human-readable Excel doc of Koh_476 / Koh_89 / COSMIC_83 classes
# with up to 10 example long_visual strings per (Koh_476, Koh_89) pair.

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
  library(ICAMS)
})

n_examples <- 10
in_dir  <- here::here("some_sup_tables")
pairs_path <- file.path(in_dir, "rosetta_stone_476_89_cap9.csv")
full_path  <- file.path(in_dir, "rosetta_stone_full_cap9.csv")
out_path   <- file.path(in_dir, "rosetta_stone_doc.xlsx")

pairs <- fread(pairs_path)
message("Loaded ", nrow(pairs), " (Koh_476, Koh_89) pairs")

message("Reading ", basename(full_path), " ...")
full <- fread(full_path)

is_single_tc_class <- function(x) {
  grepl("^[ACGT]\\[(Del|Ins)\\([CT]\\):R[^]]+\\][ACGT]$", x)
}

# For single-base T/C indels, restrict examples to ones where the indel base
# on the reported strand is itself C or T (long_visual middle token like
# "<C>" or "<T>"), and keep only one such example per Koh_476 class.
single_tc <- is_single_tc_class(full$Koh_476)
drop <- single_tc & !grepl(" <[CT]>", full$long_visual)
message("Dropping ", sum(drop),
        " rows whose strand-flipped base hides the canonical C/T")
full <- full[!drop]

message("Selecting distinct long_visual per pair ...")
full_unique <- unique(full, by = c("Koh_476", "Koh_89", "long_visual"))
setorder(full_unique, Koh_476, Koh_89)
examples <- full_unique[, {
  n <- if (is_single_tc_class(.BY$Koh_476)) 1L else n_examples
  head(.SD, n)
}, by = .(Koh_476, Koh_89),
   .SDcols = c("COSMIC_83", "long_visual")]

cosmic_per_pair <- full_unique[, .(COSMIC_83 = paste(sort(unique(COSMIC_83)),
                                                    collapse = "; ")),
                               by = .(Koh_476, Koh_89)]

doc <- merge(pairs, examples[, .(Koh_476, Koh_89, long_visual)],
             by = c("Koh_476", "Koh_89"), all.x = TRUE, sort = FALSE)
doc <- merge(doc, cosmic_per_pair, by = c("Koh_476", "Koh_89"),
             all.x = TRUE, sort = FALSE)
doc[, example_n := seq_len(.N), by = .(Koh_476, Koh_89)]

ord <- ICAMS::catalog.row.order$ID476
doc[, .row_ord := match(Koh_476, ord)]
not_in_order <- unique(doc$Koh_476[is.na(doc$.row_ord)])
if (length(not_in_order) > 0) {
  warning("Koh_476 values not in ICAMS::catalog.row.order$ID476: ",
          paste(not_in_order, collapse = ", "))
}
setorder(doc, .row_ord, Koh_89, example_n, na.last = TRUE)
doc[, .row_ord := NULL]

doc <- doc[, .(Koh_476, Koh_89, COSMIC_83, example_n, long_visual)]

# For single-base T/C indels, the flanking sequence is uninformative — keep
# only the immediate neighbouring base on each side. Then collapse the
# remaining spaces in every long_visual value.
is_single_tc <- grepl("^[ACGT]\\[(Del|Ins)\\([CT]\\):R[^]]+\\][ACGT]$",
                      doc$Koh_476)
shrink_single_tc <- function(x) {
  # x looks like "<pre> <token> <post>"; keep last char of pre, token,
  # first char of post.
  m <- regmatches(x, regexec("^(.*) (\\S+) (.*)$", x))
  vapply(m, function(p) {
    if (length(p) != 4L) return(NA_character_)
    pre <- p[2]; tok <- p[3]; post <- p[4]
    paste(substr(pre, nchar(pre), nchar(pre)), tok, substr(post, 1, 1))
  }, character(1))
}
doc$long_visual[is_single_tc] <- shrink_single_tc(doc$long_visual[is_single_tc])
doc$long_visual <- gsub(" ", "", doc$long_visual, fixed = TRUE)

n_pairs <- uniqueN(doc, by = c("Koh_476", "Koh_89"))
counts <- doc[, .N, by = .(Koh_476, Koh_89)]
n_thin <- counts[N < n_examples, .N]
message("Pairs written: ", n_pairs)
message("Pairs with < ", n_examples, " examples: ", n_thin)

# ---- Excel ----
wb <- createWorkbook()
addWorksheet(wb, "rosetta")
writeData(wb, "rosetta", doc, headerStyle = createStyle(
  textDecoration = "bold", halign = "left", border = "bottom"))

freezePane(wb, "rosetta", firstRow = TRUE)

mono <- createStyle(fontName = "Courier New")
addStyle(wb, "rosetta", mono,
         rows = seq_len(nrow(doc)) + 1,
         cols = which(names(doc) == "long_visual"),
         gridExpand = TRUE, stack = TRUE)

# Vertically merge Koh_476, Koh_89, COSMIC_83 across each pair's block.
block_cols <- match(c("Koh_476", "Koh_89", "COSMIC_83"), names(doc))
key <- paste(doc$Koh_476, doc$Koh_89, sep = "\r")
runs <- rle(key)
ends <- cumsum(runs$lengths)
starts <- ends - runs$lengths + 1L
for (i in seq_along(runs$lengths)) {
  if (runs$lengths[i] < 2) next
  rng <- (starts[i]:ends[i]) + 1L  # +1 for header row
  for (col in block_cols) mergeCells(wb, "rosetta", cols = col, rows = rng)
}

merged_style <- createStyle(valign = "top", halign = "left",
                            border = c("top", "bottom"),
                            borderColour = "grey70")
addStyle(wb, "rosetta", merged_style,
         rows = seq_len(nrow(doc)) + 1,
         cols = block_cols, gridExpand = TRUE, stack = TRUE)

setColWidths(wb, "rosetta",
             cols = seq_along(doc),
             widths = c(22, 28, 18, 8, 90))

saveWorkbook(wb, out_path, overwrite = TRUE)
message("Wrote ", out_path)
