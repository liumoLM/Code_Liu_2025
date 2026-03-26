# Plan: Update and re-run make_connection_table.R

## Context
The connection table is out of date — its exemplars don't match the current assignment file. Need to update the script to use `read_finalized()` from vhelpers.R, then re-run it.

## Changes to `Manuscript_data/Mo_CAP9_analysis/make_connection_table.R`

1. Source vhelpers.R (which provides `finalized_dir`, `finalized_files`, `read_finalized()`)
2. Remove local `finalized_dir` definition (line 6-8)
3. Replace `read.table()` calls for signatures (lines 11-22) with `read_finalized()`:
   - `sigs_83` → `as.matrix(read_finalized("83_signatures"))`
   - `sigs_89` → `as.matrix(read_finalized("89_signatures"))`
   - `sigs_476` → `as.matrix(read_finalized("476_signatures"))`
4. Replace hardcoded spectra paths (lines 25-27) with `finalized_files` entries:
   - `spectra_83_path` → `finalized_files[["83_spectra"]]`
   - `spectra_89_path` → `finalized_files[["89_spectra"]]`
   - `spectra_476_path` → `finalized_files[["476_spectra"]]`
5. Replace `out_path` (line 113) with `finalized_files[["connection_table"]]`

Then run the script with `Rscript Manuscript_data/Mo_CAP9_analysis/make_connection_table.R`.

## Verification
Check that the new connection_table.tsv has updated BestMatch89_1 values (e.g. InsDel1c should no longer have SP112907).
