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
