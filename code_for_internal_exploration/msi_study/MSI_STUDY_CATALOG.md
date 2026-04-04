# MSI Study Catalog

Brief descriptions of each R and Quarto file in `msi_study/`.

| File | Description |
|------|-------------|
| `msi_study.qmd` | Parameterized Quarto report that finds spectra similar to a given signature, validates them against annotated VCFs, plots poly-T/poly-C tract length distributions, and investigates specific mutation types of interest. |
| `render_msi_studies.R` | Driver script that renders `msi_study.qmd` multiple times with different parameters (signature, cosine cutoff, indel type) for many signatures including InsDel7, InsDel_J, InsDel_N, InsDel_P, and all signatures from the 89-to-83 mapping table. |
| `Rexemplars_7_and_J.R` | Finds exemplar spectra with cosine >= 0.9 for InsDel_J and InsDel7 (476-type) and saves the hit lists as CSV files (`J_hits_ge_0.0.csv`, `7_hits_ge_0.0.csv`). |
| `Rget_repeat_spectra.R` | For samples similar to InsDel7 (PCAWG) and InsDel_J (FMH), tallies repeat-count distributions for mono-nucleotide indels and Del2:U1:R(5,9), producing summary Excel files and detail CSVs. |
| `check_pcawg_vs_fmh_repeat_len.R` | Analyzes repeat-count distributions for the top N samples (by file size) from PCAWG and FMH datasets, producing per-pattern histograms and Del2:U1:R(5,9) detail CSVs. |
| `check_pcawg_vs_fmh_colon_repeat_len.R` | Same repeat-count analysis as above but restricted to colon/uterus/prostate cancer samples with high indel counts (>= 14000), run on PCAWG, FMH, and PCAWG graylist datasets. |
| `join_pcawg_fmh_Del2.R` | Joins PCAWG and FMH Del2:U1:R(5,9) detail CSVs, computing proportions and PCAWG/FMH proportion ratios, outputting Excel files (by short_visual+R and by R only). |
| `translate_fmh_to_pcawg.R` | Scales the top 2000 FMH 476-type spectra to approximate PCAWG-like repeat distributions by applying empirical correction factors to R(9,) rows and Del2:U1:R(5,9). |
| `Rconvert_fmh_j_to_7.R` | Similar to `translate_fmh_to_pcawg.R` but specifically for FMH samples similar to InsDel_J: scales their 476-type spectra to approximate PCAWG InsDel7-like repeat distributions. |
| `msi_study_fmh_insdel_j_converted_to_pcawg.qmd` | Quarto report that takes the FMH InsDel_J spectra (after scaling to PCAWG-like repeats) and checks how similar they are to InsDel7, with exemplar plots. |
| `clip_study/clip_study.qmd` | Like `msi_study.qmd` but uses clipped spectra (R <= 9) instead of full spectra, filtering VCFs to R <= 9 before computing catalogs. |
| `clip_study/render_clip_studies.R` | Driver script that renders `clip_study.qmd` for many signatures across 476-type, 89-type, and 83-type classifications. |
