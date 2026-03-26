# Plan: Update vignette to use new CAP9 finalized data

## Context

The finalized CAP9 analysis produced new signatures, spectra, and a new connection table in `Manuscript_data/Mo_CAP9_analysis/Finalized result/`. The vignette needs to switch to these new data files. Key structural changes: (1) the connection table now has multiple best-match exemplars per classification type instead of one curated example, (2) spectra column names no longer have cancer-type prefixes, (3) 476-type signatures now share names with 89-type (no `_476` suffix), (4) some signatures were added/removed.

## User decisions

- **Exemplar strategy**: BestMatch89_1 for 89-type views. For 83/476 views, show BOTH the type-specific best match (BestMatch83_1 / BestMatch476_1) AND BestMatch89_1's spectrum re-classified into that type.
- **Assignment file**: Keep old file; strip cancer-type prefix from its column names to match new plain sample IDs. Skip decomposition for signatures not found in assignments.
- **COSMIC/Jin/Koh comparisons**: Keep them, using reference files from `Manuscript_data/`.
- **No file moves**: All data stays where it is.

## File mapping (old → new)

| Variable | Old path | New path |
|----------|----------|----------|
| connect_89_to_83 | `Manuscript_data/89type_to_83type_connection.tsv` | `Manuscript_data/Mo_CAP9_analysis/Finalized result/connection_table.tsv` |
| type83_spectra | `Manuscript_data/Liu_et_al_83_type_spectra.tsv` | `Finalized result/liu_et_al_83_spectra.tsv` |
| type83_sigs | `Manuscript_data/Liu_et_al_final_83_type_signatures.tsv` | `Finalized result/liu_et_al_83_signatures.tsv` |
| type89_sigs | `Manuscript_data/Liu_et_al_final_89_type_signatures.tsv` | `Finalized result/liu_et_al_89_signatures.tsv` |
| type89_spectra | `Manuscript_data/Liu_et_al_89_type_spectra.tsv` | `Finalized result/liu_et_al_89_spectra.tsv` |
| type476_sigs | `Manuscript_data/Liu_et_al_final_476_type_signatures.tsv` | `Finalized result/liu_et_al_476_signatures.tsv` |
| ID476_spectra | `Manuscript_data/Liu_et_al_476_type_spectra.tsv` | `Finalized result/liu_et_al_476_spectra.tsv` |
| cosmic_sigs | `Manuscript_data/COSMIC_v3.5_ID_GRCh37_signatures.tsv` | **unchanged** |
| jin_sigs | `Manuscript_data/jin_2024_sup_tab_1_signatures.tsv` | **unchanged** |
| koh_sigs | `Manuscript_data/Koh_signatures.tsv` | **unchanged** |
| assignment | `Manuscript_data/Liu_et_al_89_type_signature_assignments.tsv` | **unchanged** (strip cancer-type prefix from colnames) |

## Changes by file

### 1. `vignette/vignette.qmd`

**a. Add finalized data directory** (after line 61)
```r
finalized_dir <- file.path(data_dir, "Mo_CAP9_analysis", "Finalized result")
```
Keep `data_dir` pointing to `Manuscript_data/` for COSMIC/Jin/Koh/assignment files.

**b. Connection table (lines 109-122)**
- Read from `finalized_dir/connection_table.tsv`
- No column renaming needed (keep original names: `InDel89`, `InDel83`, `BestMatch83_1`, etc.)
- OR rename to: `ID89_signature`, `ID83_signature`, `exemplar_83`, ..., `exemplar_89`, ..., `exemplar_476`, ...

**c. Spectra/signature reads (lines 129-354)**
- All 6 data files: read from `finalized_dir` with new filenames
- New files have **no `rowname` column** — mutation types are already row names. Use `row.names = 1` since the first column (unnamed) contains the row names.
- Remove `data.table::fread` workarounds where they were used to handle the old format; or keep fread and adjust.

**d. Assignment file (lines 327-334)**
- Keep reading from `data_dir` (old location)
- After reading, strip cancer-type prefix from column names:
  ```r
  colnames(ID89.mSigAct.assignment) <- sub("^.*::", "", colnames(ID89.mSigAct.assignment))
  ```

**e. `compute_sig_data` call (lines 408-433)**
- Pass three exemplar IDs per signature from connection table:
  ```r
  exemplar_89 = connect_89_to_83$BestMatch89_1[i],
  exemplar_83 = connect_89_to_83$BestMatch83_1[i],
  exemplar_476 = connect_89_to_83$BestMatch476_1[i]
  ```

### 2. `vignette/vhelpers.R`

**a. `compute_sig_data()` (lines 87-279)**
- Change signature: replace `exemplar_id` with `exemplar_89`, `exemplar_83`, `exemplar_476`
- Store all three in result list (for use by onesig.qmd)
- **cosine89**: sig89 vs `ID89_catalogs[, exemplar_89]`
- **cosine83**: sig83 vs `ID83_catalogs[, exemplar_83]`; ALSO compute sig83 vs `ID83_catalogs[, exemplar_89]` → store as `cosine83_linking`
- **cosine476**: sig476 vs `ID476_catalogs[, exemplar_476]`; ALSO compute sig476 vs `ID476_catalogs[, exemplar_89]` → store as `cosine476_linking`
- **Decomposition**: Use `exemplar_89` for the residual/partial spectrum computation. Check if `exemplar_89` exists in assignment matrix columns; if not, skip decomposition gracefully.
- Update `is_polyT_removed` if signature names changed
- The `InsDel_N → InsDel_J` mapping: check if still needed with new signature set (InsDel_N_alpha/beta in new data)

**b. `generate_plots_to_files()` (lines 324-800)**
- Add new plot paths: `id83_catalog_83match`, `id476_catalog_476match` (the type-specific best match exemplar plots)
- For 83-type section: generate plot for `exemplar_83`'s 83-type spectrum (the native best match) AND `exemplar_89`'s 83-type spectrum
- For 476-type section: generate plot for `exemplar_476`'s 476-type spectrum (the native best match) AND `exemplar_89`'s 476-type spectrum
- For 89-type section: use `exemplar_89` (current behavior, just different variable name)

**c. `reconstruct_plot_paths()` (lines 999-1183)**
- Add reconstruction logic for new plot types (83-match catalog, 476-match catalog)

**d. `check_plot_cache()` / `save_plot_cache()` (lines 888-989)**
- Update file name list to match new filenames from finalized_dir

### 3. `vignette/onesig.qmd`

**a. 476-type section (lines 112-167)**
- After existing 476 catalog plot (BestMatch89_1's 476-type spectrum), add BestMatch476_1's 476-type spectrum with cosine similarity text

**b. 83-type section (lines 291-335)**
- After existing 83 catalog plot (BestMatch89_1's 83-type spectrum), add BestMatch83_1's 83-type spectrum with cosine similarity text
- Update text that refers to `sig_data$exemplar_id` → `sig_data$exemplar_89` (or similar)

**c. Update references to `exemplar_id`**
- Lines 133, 157-159, 311-312: change `sig_data$exemplar_id` to appropriate exemplar variable

### 4. `vignette/onesig_standalone.qmd`

- Update `data_dir` to handle both directories (or add `finalized_dir`)
- Minor: if vhelpers.R is also sourced here, changes propagate automatically

## Signature set changes

Old (48 rows in connection table): includes InsDel1d, InsDel5a, InsDel5b, InsDel_K_alpha, InsDel_K_beta, InsDel_O, InsDel_P
New (46 rows): includes InsDel5c, InsDel_K, InsDel_N_alpha, InsDel_N_beta; removes InsDel1d, InsDel_K_alpha, InsDel_K_beta, InsDel_O, InsDel_P

The `is_insdel15_16` check stays (InsDel15/16 still present).
The `is_polyT_removed` check: currently checks for `C_ID7, ID_J, C_ID10, ID_N, ID_O`. Since ID_O is gone from new data, remove from list. ID_N may need update if InsDel_N_alpha/beta map differently.

## Verification

1. Run `quarto render vignette.qmd` and check for errors
2. Verify all ~46 signature sections render with correct plots
3. Check that 83-type and 476-type sections show two exemplar spectra each
4. Verify decomposition plots appear for signatures that exist in old assignment file
5. Verify COSMIC/Jin/Koh comparison plots still render correctly
6. Spot-check a few cosine similarities
