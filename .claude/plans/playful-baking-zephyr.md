# Plan: Move find_many_similar.R and create find_similar_to_insdel_N_beta.R

## Context
Moving `find_many_similar.R` from `old_code/` to `code_for_internal_exploration/` and writing a script that uses it to find spectra similar to InsDel_N_beta across both 89-type and 476-type classifications.

## Steps

### 1. Move file
- `git mv old_code/find_many_similar.R code_for_internal_exploration/find_many_similar.R`

### 2. Create output directory
- `mkdir -p code_for_internal_exploration/similar_to_indesdel_N_beta/`

### 3. Create `code_for_internal_exploration/similar_to_indesdel_N_beta/find_similar_to_insdel_N_beta.R`

Short script that:
- Sources `../find_many_similar.R`
- Runs `find_many_similar()` twice with `do_plot = FALSE`:
  1. **476-type**: sig_path = `Manuscript_data/finalized_cap9/liu_et_al_476_signatures.tsv`, spectra_path = `Manuscript_data/finalized_cap9/liu_et_al_476_spectra.tsv`, sig_col = `"InsDel_N_beta"`
  2. **89-type**: sig_path = `Manuscript_data/finalized_cap9/liu_et_al_89_signatures.tsv`, spectra_path = `Manuscript_data/finalized_cap9/liu_et_al_89_spectra.tsv`, sig_col = `"InsDel_N_beta"`
- Parameters: `cosine_cutoff = 0.9`, `num_exemplars = 300`
- Saves results (the `above_cutoff` data frames) as TSVs in the same directory

## Files modified
- `old_code/find_many_similar.R` → `code_for_internal_exploration/find_many_similar.R` (move)
- `code_for_internal_exploration/similar_to_indesdel_N_beta/find_similar_to_insdel_N_beta.R` (new)

## Verification
- Run the script from the project root and confirm the two output TSVs are created
