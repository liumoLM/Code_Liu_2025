# Plan: Format Jin matches like Koh in overview tables

## Context
Jin signatures in `prot_table_1.csv`/`.xlsx` and the HTML overview table should be treated like Koh signatures: don't show matches below 0.9 (already done via `Jin_min_cosine <- 0.9`), and when the same Jin signature matches multiple of our signatures, make only the best match (highest cosine) dark and gray out the rest.

The Koh duplicate-handling pattern already exists in two files and should be replicated for Jin.

## Changes

### 1. `vignette/table_1_as_dt.R` — HTML overview table

After the existing Koh duplicate-handling block (lines 31-59), add an identical block for Jin:
- Operate on `best_match_jin` and `cosine_v_jin` columns
- Find duplicates, mark best with `*`, collect non-best rows
- Wrap non-best rows in gray `<span>` tags

Also round `cosine_v_jin` before the HTML wrapping (like `cos_v_koh` on line 29).

### 2. `vignette/prot_table_1_to_excel.R` — Excel export

After the existing Koh duplicate-handling block (lines 25-51), add an identical block for Jin:
- Operate on `best_match_jin` and `cosine_v_jin` columns
- Track `jin_rows_to_gray`

After the existing Koh gray-styling block (lines 103-132), add gray styling for `jin_rows_to_gray` on the Jin columns.

### No changes needed to `vignette.qmd`
The 0.9 threshold already filters `jin_matches` upstream, so `best_match_jin` will be NA for sigs without a 0.9+ match.

## Verification
- Run `Rscript prot_table_1_to_excel.R` and inspect the Excel file for gray Jin duplicates
- Render the vignette and check the HTML overview table shows gray Jin duplicates
