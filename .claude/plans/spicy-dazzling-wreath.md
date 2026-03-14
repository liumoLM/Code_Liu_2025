# Plan: Parameterize cluster_cap9_koh89 by mutation_type (89 or 476)

## Context
The clustering QMD and its render script currently only handle 89-type signatures. We need to add a `mutation_type` parameter (default 89) so the same pipeline works for 476-type data.

## Files to modify
1. `Manuscript_data/Mo_CAP9_analysis/render_cluster.R`
2. `Manuscript_data/Mo_CAP9_analysis/cluster_cap9_koh89.qmd`

## Input data files

**When mutation_type = 89 (current behavior):**
- Signatures: `Manuscript_data/Mo_CAP9_analysis/Signatures/Koh89/CAP9.*.txt`
- Catalogs: `Manuscript_data/Mo_CAP9_analysis/Catalogs/CAP9.PCAWG.Koh89.catalog.txt` and `CAP9.Hartwig.Koh89.catalog.txt`
- Liu reference: `Manuscript_data/Liu_et_al_final_89_type_signatures.tsv`
- Koh reference: `Manuscript_data/Koh_signatures.tsv`

**When mutation_type = 476:**
- Signatures: `Manuscript_data/Mo_CAP9_analysis/Signatures/Koh476/CAP9.*.txt`
- Catalogs: `Manuscript_data/Mo_CAP9_analysis/Catalogs/CAP9.PCAWG.Koh476.catalog.txt` and `CAP9.Hartwig.Koh476.catalog.txt`
- Liu reference: `Manuscript_data/Liu_et_al_final_476_type_signatures.tsv`
- Koh reference: **none** (no 476-type Koh reference exists)

## Changes to render_cluster.R

Add `--mutation-type` argument (numeric, default 89). Pass it to `quarto_render()` as `execute_params`. Update output filename to use `mutation_type` instead of hardcoded `89`:
- `cluster_cap9_89_minsim_0.95.html` → `cluster_cap9_{mutation_type}_minsim_{sim}.html`

## Changes to cluster_cap9_koh89.qmd

### 1. YAML params (line 7)
Add `mutation_type: 89`.

### 2. `params` chunk (lines 40–45)
- Read `params$mutation_type` into `mutation_type`
- Set `sig_dir` dynamically: `Signatures/Koh89` or `Signatures/Koh476`
- Set `koh_label` to `"Koh89"` or `"Koh476"` for use in glob patterns and file names

### 3. Title (line 2)
Include `mutation_type` in the title dynamically.

### 4. `load-sigs` chunk (lines 49–82)
- Change hardcoded `nrows == 89` filter to `nrows == mutation_type`
- Change `*Koh89*.txt` catalog glob to use `koh_label`
- Change `ICAMS::catalog.row.order$ID89` to `ICAMS::catalog.row.order[[paste0("ID", mutation_type)]]`

### 5. `medoid-ref-load` chunk (lines 229–257)
- Liu reference: switch file path based on `mutation_type` (`Liu_et_al_final_89_type_signatures.tsv` vs `Liu_et_al_final_476_type_signatures.tsv`)
- Koh reference: only load `Koh_signatures.tsv` when `mutation_type == 89`; skip for 476
- Adjust `ref_combined_med`, `ref_combined_all`, `ref_source_med`, `ref_source_all` to conditionally include/exclude Koh

### 6. `load-catalogs` chunk (lines 540–558)
- Switch catalog filenames: `CAP9.PCAWG.Koh89.catalog.txt` → `CAP9.PCAWG.{koh_label}.catalog.txt` (same for Hartwig)

### 7. `medoid-plots` chunk (lines 560–617)
- Replace `plot_89()` calls with conditional: `plot_89()` for 89, `plot_476()` for 476

### 8. `save-outputs` chunk (lines 621–656)
- Change prefix from hardcoded `"CAP9_Koh89_clustered"` to `paste0("CAP9_Koh", mutation_type, "_clustered")`

## Verification
- Render 89-type: `Rscript Manuscript_data/Mo_CAP9_analysis/render_cluster.R` (default, should work as before)
- Render 476-type: `Rscript Manuscript_data/Mo_CAP9_analysis/render_cluster.R --mutation-type 476`
- Check both produce correct HTML output files
