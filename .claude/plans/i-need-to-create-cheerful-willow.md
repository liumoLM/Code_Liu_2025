# Human-readable Koh_476 / Koh_89 / long_visual documentation

## Context

`some_sup_tables/rosetta_stone_476_89_cap9.csv` (473 unique Koh_476–Koh_89 pairs)
and `some_sup_tables/rosetta_stone_full_cap9.csv` (~16 M raw observations with
columns `Koh_476, Koh_89, COSMIC_83, long_visual`) together let a reader see
*examples* of the genomic context that each Koh_476 / Koh_89 class summarises.
Right now those examples are buried in a 1.6 GB CSV. We want a single,
flippable artefact that, for each Koh_476 / Koh_89 pair, shows up to 10
representative `long_visual` strings — useful as supplementary documentation
for the paper and for anyone trying to learn the indel naming convention.

## Approach

Write a new script `some_sup_tables/build_rosetta_doc.R` that:

1. Loads the 473-pair table `rosetta_stone_476_89_cap9.csv`.
2. Streams `rosetta_stone_full_cap9.csv` with `data.table::fread`. Groups by
   `(Koh_476, Koh_89)` and, per group, takes up to **10 distinct**
   `long_visual` values (preserving file order — these are already a mix of
   FMH + PCAWG samples). Distinct rather than random so the table shows
   variety, not duplicates.
3. Joins the per-pair long_visual list back to the 473-pair table. Also
   carries `COSMIC_83` (collapsed to the unique set per pair — usually one
   value, sometimes a small set). Including COSMIC_83 makes the doc more
   useful and costs nothing.
4. Orders rows by `match(Koh_476, ICAMS::catalog.row.order$ID476)`. Any pair
   not in that vector (none expected, but defensive) sorts to the end.
5. Expands to one row per `long_visual` example (≤10 rows per pair) and
   writes an `.xlsx` via **openxlsx** (the package already used in this repo
   — `vignette/`, `msi_study/`).

## Output format (recommended)

`some_sup_tables/rosetta_stone_doc.xlsx` — one sheet with columns

| Koh_476 | Koh_89 | COSMIC_83 | example_# | long_visual |

Within each `(Koh_476, Koh_89)` block (≤10 rows):
- The `Koh_476`, `Koh_89`, and `COSMIC_83` cells are **merged vertically**
  (`openxlsx::mergeCells`) so each pair occupies one visual row spanning its
  examples.
- A thin border between blocks (`openxlsx::addStyle`) for scannability.
- Header row bold, frozen pane on row 1 (`freezePane`).
- Column widths: Koh_476 / Koh_89 / COSMIC_83 sized to fit, `long_visual`
  ~80 chars wide (monospace font for that column so the `<…>` flank
  visualisation stays aligned).

This is what the user explicitly suggested ("merged to one row each for
Koh_476 and Koh_89") and matches the established repo convention.

I considered two alternatives but am not recommending them:
- One row per pair with `long_visual` examples newline-joined in a single
  cell. Compact, but Excel's wrap-text rendering breaks the monospace
  alignment of the visual flanks.
- Plain CSV with `;`-joined examples. Loses the visual benefit entirely.

## Critical files

- **New**: `some_sup_tables/build_rosetta_doc.R`
- **New output**: `some_sup_tables/rosetta_stone_doc.xlsx`
- **Reads**: `some_sup_tables/rosetta_stone_476_89_cap9.csv`,
  `some_sup_tables/rosetta_stone_full_cap9.csv`
- **Reuses**: `ICAMS::catalog.row.order$ID476` (character(476) indel ordering),
  `openxlsx` (already a project dependency), `data.table::fread` (used elsewhere
  in this folder, e.g. `cap9_filter.R`).

## Verification

```
Rscript some_sup_tables/build_rosetta_doc.R
xdg-open some_sup_tables/rosetta_stone_doc.xlsx
```

Sanity checks the script should print:
- Number of (Koh_476, Koh_89) pairs written (should equal 473).
- Number of pairs with fewer than 10 distinct `long_visual` examples (rare
  classes — flag them so we know the doc is thinner there).
- Any Koh_476 in the rosetta that is NOT in `ICAMS::catalog.row.order$ID476`
  (should be 0; warn if not).

---

## Iteration 2 — formatting polish + rich-text long_visual

Stepwise refinements requested after the first build:

1. **20 examples per Koh_476** (was 10). For single-base T/C indels we keep
   the previous one-row rule. So `n_examples` becomes 20 only for the
   non-single-T/C classes.
2. **Vertical-center every row** including the merged Koh_476 / Koh_89 /
   COSMIC_83 cells (`valign = "center"`, replacing the current `"top"`).
3. **Wrap text in column C** (`COSMIC_83`). Some pairs concatenate multiple
   COSMIC types with `"; "` and need wrapping to stay readable. Use
   `createStyle(wrapText = TRUE, valign = "center")`.
4. **Asterisk on ambiguous Koh_476**. A few Koh_476 classes (e.g.
   `Del(6,):U(4,):R(2,9)`) appear in two rows because they map to two
   different Koh_89 values. Compute
   `multi <- doc[, uniqueN(Koh_89), by = Koh_476][V1 > 1]$Koh_476` and append
   `"*"` to the Koh_476 cell wherever it occurs. List the asterisked classes
   in the script's stdout.
5. **Rich-text colouring inside `long_visual`**:
   - characters between `<` and `>` → **bold + blue**
   - characters between `[` and `]` → **blue** (not bold)
   - characters between `{` and `}` → **underlined** (default colour)
   - the bracket characters themselves keep the default style.

### Why this needs an extra package

`openxlsx` (currently used) does **not** support cell-level rich text.
`openxlsx2` does, via `fmt_txt()` chunks plus `wb_add_data`. `openxlsx2`
is not currently installed (`requireNamespace("openxlsx2") == FALSE`).

The cheapest path is:

- Install `openxlsx2` (`install.packages("openxlsx2")`). Pure-R + Rcpp
  dependency, no system libs needed.
- Migrate `build_rosetta_doc.R` from `openxlsx` to `openxlsx2`. The rest of
  this repo's `openxlsx` usage stays untouched.

If installing isn't acceptable, the fallback is to drop the within-cell
colour requirement (still applying bold-blue / blue / underline to the
*entire* `long_visual` cell when it contains those bracket types, which is
much less informative).

### Implementation sketch

```r
library(openxlsx2)

# helper: split "X<abc>Y[def]Z{ghi}W" into a list of fmt_txt() chunks with
# the right styling applied to the inner characters of each bracket type.
rich_long_visual <- function(s) {
  # tokenise: text | <...> | [...] | {...}
  pieces <- regmatches(
    s,
    gregexpr("<[^>]*>|\\[[^]]*\\]|\\{[^}]*\\}|[^<\\[\\{]+", s, perl = TRUE)
  )[[1]]
  out <- lapply(pieces, function(p) {
    if (startsWith(p, "<")) {
      c(fmt_txt("<"),
        fmt_txt(substr(p, 2, nchar(p) - 1), bold = TRUE, color = wb_color("blue")),
        fmt_txt(">"))
    } else if (startsWith(p, "[")) {
      c(fmt_txt("["),
        fmt_txt(substr(p, 2, nchar(p) - 1), color = wb_color("blue")),
        fmt_txt("]"))
    } else if (startsWith(p, "{")) {
      c(fmt_txt("{"),
        fmt_txt(substr(p, 2, nchar(p) - 1), underline = TRUE),
        fmt_txt("}"))
    } else {
      fmt_txt(p)
    }
  })
  Reduce(`+`, unlist(out, recursive = FALSE))
}

wb <- wb_workbook()$add_worksheet("rosetta")
wb$add_data(x = doc[, .(Koh_476, Koh_89, COSMIC_83, example_n)])
# write long_visual cell-by-cell as rich text
for (i in seq_len(nrow(doc))) {
  wb$add_data(x = rich_long_visual(doc$long_visual[i]),
              dims = wb_dims(rows = i + 1, cols = 5))
}
# styles: vertical center everywhere, wrapText for column C, monospace
# Courier on column 5, freeze pane, merge Koh_476/Koh_89/COSMIC_83 across
# pair blocks (openxlsx2's wb_merge_cells), bold header, column widths.
wb_save(wb, out_path)
```

(Real implementation will iterate cells in a single `add_data` batch where
possible and replicate the existing merge/border logic.)

### Updated critical files

- **Modify**: `some_sup_tables/build_rosetta_doc.R` (port from `openxlsx`
  → `openxlsx2`, add the formatting and asterisk logic, bump `n_examples`).
- **New runtime dependency**: `openxlsx2` (CRAN; add to user library via
  `install.packages("openxlsx2")`).

### Updated verification

```
Rscript -e 'install.packages("openxlsx2", repos = "https://cloud.r-project.org")'
Rscript some_sup_tables/build_rosetta_doc.R
xdg-open some_sup_tables/rosetta_stone_doc.xlsx
```

Visual checks in the spreadsheet:
- Koh_476 / Koh_89 / COSMIC_83 cells centered both horizontally and
  vertically.
- COSMIC_83 cells wrap when they hold multiple values.
- Koh_476 values that map to >1 Koh_89 end with `*` (script also lists them
  to stdout).
- In `long_visual`: `<…>` shows bold blue inside, `[…]` plain blue inside,
  `{…}` underlined inside.
- 20 example rows for typical pairs; still 1 for single-base T/C classes.
