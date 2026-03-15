# Plan: Convert 476-type medoid signatures to 89-type

## Context
User needs a function that reads the 476-type clustered medoid signatures and converts them to 89-type using `t476_to_89()`, saving the result in the same directory.

## Implementation

Create `code/collapse_for_dendro.R` with a function `collapse_for_dendro()` that:

1. Sources `code/collapse_476_to_89.R` to get `t476_to_89()`
2. Reads `Manuscript_data/Mo_CAP9_analysis/clustering_results/CAP9_Koh476_clustered_medoid_signatures.tsv`
3. Calls `t476_to_89()` on the data frame
4. Strips the `_converted` suffix from column names
5. Writes the result to `Manuscript_data/Mo_CAP9_analysis/clustering_results/CAP9_476_converted_to_89.tsv`

### Key details
- **Input**: TSV with 476 mutation type row names, 70 signature columns
- **`t476_to_89()`** appends `_converted` suffix to column names — must strip it
- **Output**: 89 rows × 70 columns, original column names, written with row names

### Files
- **Create**: `code/collapse_for_dendro.R`
- **Read (reuse)**: `code/collapse_476_to_89.R` (sourced for `t476_to_89`)

## Verification
- Run `source("code/collapse_for_dendro.R"); collapse_for_dendro()`
- Check output has 89 rows, 70 columns
- Verify row names match `ICAMS::catalog.row.order$ID89`
- Verify column names match original input (no `_converted` suffix)
