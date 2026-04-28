#!/usr/bin/env Rscript
# Build a human-readable Excel doc of Koh_476 / Koh_89 / COSMIC_83 classes
# with example long_visual strings per (Koh_476, Koh_89) pair.

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx2)
  library(ICAMS)
})

n_examples <- 20 # for non single-base T/C indels
n_examples_singletc <- 1
in_dir <- here::here("some_sup_tables")
pairs_path <- file.path(in_dir, "rosetta_stone_476_89_cap9.csv")
full_path <- file.path(in_dir, "rosetta_stone_full_cap9.csv")
out_path <- file.path(in_dir, "sup_table_rosetta_stone.xlsx")

is_single_tc_class <- function(x) {
  grepl("^[ACGT]\\[(Del|Ins)\\([CT]\\):R[^]]+\\][ACGT]$", x)
}

pairs <- fread(pairs_path)
message("Loaded ", nrow(pairs), " (Koh_476, Koh_89) pairs")

message("Reading ", basename(full_path), " ...")
full <- fread(full_path)

# For single-base T/C indels, restrict examples to ones where the indel base
# on the reported strand is itself C or T (long_visual middle token like
# "<C>" or "<T>"), and keep only one such example per Koh_476 class.
single_tc <- is_single_tc_class(full$Koh_476)
drop <- single_tc & !grepl(" <[CT]>", full$long_visual)
message(
  "Dropping ",
  sum(drop),
  " rows whose strand-flipped base hides the canonical C/T"
)
full <- full[!drop]

message("Selecting distinct long_visual per pair ...")
full_unique <- unique(full, by = c("Koh_476", "Koh_89", "long_visual"))
setorder(full_unique, Koh_476, Koh_89)
examples <- full_unique[,
  {
    n <- if (is_single_tc_class(.BY$Koh_476)) {
      n_examples_singletc
    } else {
      n_examples
    }
    head(.SD, n)
  },
  by = .(Koh_476, Koh_89),
  .SDcols = c("COSMIC_83", "long_visual")
]

cosmic_per_pair <- full_unique[,
  .(COSMIC_83 = paste(sort(unique(COSMIC_83)), collapse = "; ")),
  by = .(Koh_476, Koh_89)
]

doc <- merge(
  pairs,
  examples[, .(Koh_476, Koh_89, long_visual)],
  by = c("Koh_476", "Koh_89"),
  all.x = TRUE,
  sort = FALSE
)
doc <- merge(
  doc,
  cosmic_per_pair,
  by = c("Koh_476", "Koh_89"),
  all.x = TRUE,
  sort = FALSE
)
doc[, example_n := seq_len(.N), by = .(Koh_476, Koh_89)]

# Order by ICAMS ID476 row order.
ord <- ICAMS::catalog.row.order$ID476
doc[, .row_ord := match(Koh_476, ord)]
not_in_order <- unique(doc$Koh_476[is.na(doc$.row_ord)])
if (length(not_in_order) > 0) {
  warning(
    "Koh_476 values not in ICAMS::catalog.row.order$ID476: ",
    paste(not_in_order, collapse = ", ")
  )
}
setorder(doc, .row_ord, Koh_89, example_n, na.last = TRUE)
doc[, .row_ord := NULL]

# For single-base T/C indels, shrink long_visual flanks to the immediate
# neighbours; then collapse spaces in every long_visual.
is_single_tc_doc <- is_single_tc_class(doc$Koh_476)
shrink_single_tc <- function(x) {
  m <- regmatches(x, regexec("^(.*) (\\S+) (.*)$", x))
  vapply(
    m,
    function(p) {
      if (length(p) != 4L) {
        return(NA_character_)
      }
      pre <- p[2]
      tok <- p[3]
      post <- p[4]
      paste(substr(pre, nchar(pre), nchar(pre)), tok, substr(post, 1, 1))
    },
    character(1)
  )
}
doc$long_visual[is_single_tc_doc] <- shrink_single_tc(doc$long_visual[
  is_single_tc_doc
])
doc$long_visual <- gsub(" ", "", doc$long_visual, fixed = TRUE)

# Asterisk Koh_476 values that map to more than one Koh_89.
multi_koh476 <- doc[, uniqueN(Koh_89), by = Koh_476][V1 > 1, Koh_476]
if (length(multi_koh476) > 0) {
  message(
    "Koh_476 classes mapping to multiple Koh_89 (marked with *): ",
    paste(multi_koh476, collapse = ", ")
  )
  doc[Koh_476 %in% multi_koh476, Koh_476 := paste0(Koh_476, "*")]
}

doc <- doc[, .(Koh_476, Koh_89, COSMIC_83, example_n, long_visual)]

n_pairs <- uniqueN(doc, by = c("Koh_476", "Koh_89"))
counts <- doc[, .N, by = .(Koh_476, Koh_89)]
n_thin <- counts[N < n_examples, .N]
message("Pairs written: ", n_pairs)
message("Pairs with < ", n_examples, " examples: ", n_thin)

# ---- Excel via openxlsx2 ----
wb <- wb_workbook()$add_worksheet("rosetta")

# Write the non-rich columns first.
plain <- copy(doc)
plain[, long_visual := NA_character_]
wb$add_data(x = plain, na.strings = "")

# Write long_visual as rich text per cell.
cyan <- wb_color(hex = "FFFF1493") # bright pink (DeepPink)
mono_font <- "Liberation Mono"

# Recursively style: top level recognises <...>, [...], {...}; inside <...>
# any nested {...} is additionally underlined.
rich_long_visual <- function(s) {
  if (is.na(s) || !nzchar(s)) {
    return(fmt_txt(""))
  }
  txt <- function(x, ...) fmt_txt(x, font = mono_font, ...)
  pieces <- regmatches(
    s,
    gregexpr("<[^>]*>|\\[[^]]*\\]|\\{[^}]*\\}|[^<\\[\\{]+", s, perl = TRUE)
  )[[1]]
  parts <- list()
  for (p in pieces) {
    if (startsWith(p, "<")) {
      inner <- substr(p, 2L, nchar(p) - 1L)
      parts <- c(parts, list(txt("<")))
      subs <- regmatches(
        inner,
        gregexpr("\\{[^}]*\\}|[^{]+", inner, perl = TRUE)
      )[[1]]
      for (sp in subs) {
        if (startsWith(sp, "{")) {
          parts <- c(
            parts,
            list(txt("{")),
            list(txt(
              substr(sp, 2L, nchar(sp) - 1L),
              bold = TRUE,
              color = cyan,
              underline = TRUE
            )),
            list(txt("}"))
          )
        } else {
          parts <- c(parts, list(txt(sp, bold = TRUE, color = cyan)))
        }
      }
      parts <- c(parts, list(txt(">")))
    } else if (startsWith(p, "[")) {
      parts <- c(
        parts,
        list(txt("[")),
        list(txt(substr(p, 2L, nchar(p) - 1L), color = cyan)),
        list(txt("]"))
      )
    } else if (startsWith(p, "{")) {
      parts <- c(
        parts,
        list(txt("{")),
        list(txt(substr(p, 2L, nchar(p) - 1L), underline = TRUE)),
        list(txt("}"))
      )
    } else {
      parts <- c(parts, list(txt(p)))
    }
  }
  Reduce(`+`, parts)
}

long_col <- which(names(doc) == "long_visual")
for (i in seq_len(nrow(doc))) {
  wb$add_data(
    x = rich_long_visual(doc$long_visual[i]),
    dims = wb_dims(rows = i + 1L, cols = long_col)
  )
}

# Header style
wb$add_cell_style(
  dims = wb_dims(rows = 1, cols = seq_along(doc)),
  horizontal = "left",
  vertical = "center"
)
wb$add_font(dims = wb_dims(rows = 1, cols = seq_along(doc)), bold = "1")

# All body cells: vertical center.
body_rows <- seq_len(nrow(doc)) + 1L
wb$add_cell_style(
  dims = wb_dims(rows = body_rows, cols = seq_along(doc)),
  vertical = "center"
)

# Wrap text in COSMIC_83 column.
wb$add_cell_style(
  dims = wb_dims(rows = body_rows, cols = which(names(doc) == "COSMIC_83")),
  wrap_text = "1",
  vertical = "center",
  horizontal = "left"
)

# Monospace long_visual cells (font also set per rich-text run, but set the
# whole-cell font too so plain segments stay aligned).
wb$add_font(dims = wb_dims(rows = body_rows, cols = long_col), name = mono_font)

wb$freeze_pane(first_row = TRUE)

# Vertically merge Koh_476 / Koh_89 / COSMIC_83 across each pair's block.
block_cols <- match(c("Koh_476", "Koh_89", "COSMIC_83"), names(doc))
key <- paste(doc$Koh_476, doc$Koh_89, sep = "\r")
runs <- rle(key)
ends <- cumsum(runs$lengths)
starts <- ends - runs$lengths + 1L
for (i in seq_along(runs$lengths)) {
  if (runs$lengths[i] < 2) {
    next
  }
  rng <- (starts[i]:ends[i]) + 1L
  for (col in block_cols) {
    wb$merge_cells(dims = wb_dims(rows = rng, cols = col))
  }
}

# Column widths.
wb$set_col_widths(cols = seq_along(doc), widths = c(22, 28, 18, 8, 90))

wb_save(wb, out_path, overwrite = TRUE)
message("Wrote ", out_path)
