# Partial credit of signatures to Del(T) from long T-homopolymers

## Context
We want to quantify each signature's contribution to deletions of T from stretches of 6, 7, 8, or 9 T's for every tumor. This uses "partial credit": for a given signature assigned N mutations to a tumor, the credit for a specific 476-type channel is `sig_profile[channel] * N`. We sum these partial credits across the 36 target Del(T) channels (R6, R7, R8, R(9,) × 9 flanking-base combos) to get one number per signature per tumor.

## Input files
- `Manuscript_data/finalized_cap9/liu_et_al_89_assignment.tsv` — 42 signatures × ~6975 samples (mutation counts)
- `Manuscript_data/finalized_cap9/liu_et_al_476_signatures.tsv` — 476 channels × 44 signatures (signature profiles, same sig names as assignment file)
- `Manuscript_data/sample_info.tsv` — Patient, Cancer_Type, etc.

## Target channels (36 total)
Regex on ID476 row names: `Del\(T\):R(6|7|8|\(9,\))`
- `[ACG][Del(T):R6][ACG]` — 9 channels
- `[ACG][Del(T):R7][ACG]` — 9 channels
- `[ACG][Del(T):R8][ACG]` — 9 channels
- `[ACG][Del(T):R(9,)][ACG]` — 9 channels

## New script: `code_for_internal_exploration/partial_credit_delT_long_homopolymers.R`

### Algorithm
1. Load assignment matrix (sigs × samples), 476-type signatures (channels × sigs), and sample_info
2. Identify the 36 target channels via `grep("Del\\(T\\):R(6|7|8|\\(9,\\))", rownames(sigs_476))`
3. Extract the sub-matrix of 476 signatures for just those 36 channels: `target_profile` (36 × n_sigs)
4. For each signature, sum the target profile values: `target_fraction[sig] = colSums(target_profile)[sig]`
5. Compute partial credit matrix: for each signature s and sample t, `partial_credit[s, t] = target_fraction[s] * assignment[s, t]`
   - This is a simple element-wise multiplication: `partial_credit = diag(target_fraction) %*% assignment` or vectorized `sweep(assignment, 1, target_fraction, "*")`
6. Result is a matrix (sigs × samples) of partial credits for Del(T) from long T-homopolymers
7. Transpose and join with sample_info to add Cancer_Type
8. Write output TSV

### Output
- TSV file with columns: Patient, Cancer_Type, then one column per signature containing partial credit values
- Output to `plot_output/` or `/tmp/`

## Verification
- Run the script: `Rscript code_for_internal_exploration/partial_credit_delT_long_homopolymers.R`
- Spot-check: for a sample with known high InsDel7 (MSI signature), the InsDel7 partial credit should be large
- Verify column sums make sense relative to total Del(T) counts in the 476 spectra
