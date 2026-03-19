# Rename "Finalized result" → "finalized_cap9" in all references

## Context
The directory `Manuscript_data/Mo_CAP9_analysis/Finalized result/` was renamed to `finalized_cap9/`. All code references need updating.

## Files to edit (5 files, 1 change each)

1. **`vignette/vignette.qmd`** line 62: `"Finalized result"` → `"finalized_cap9"`
2. **`vignette/onesig_standalone.qmd`** line 31: same
3. **`vignette/render_separate_pages.R`** line 107: same
4. **`code/collapse_476_to_83.R`** line 810: same
5. **`Manuscript_data/Mo_CAP9_analysis/make_connection_table.R`** line 7: same

## Verification
- `grep -r "Finalized result" --include="*.R" --include="*.qmd"` should return no hits (excluding `.claude/plans/`)
