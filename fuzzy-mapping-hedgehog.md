# Plan: Implement MSI Study Quarto Report

## Context
Create a parameterized Quarto document (`msi_study/msi_study.qmd`) that reports on MSI tumor characteristics for a given indel signature: poly-T/C tract length distributions, and detailed breakdowns of specific mutation types of interest. The spec is in `msi_study/msi_study_specification.md`.

## Parameters (YAML params block)
- `sig_to_report`: signature ID (e.g. `"InsDel_J"`)
- `types_of_interest`: list of 476-type mutation types (default: `["Del2:U1:R(5,9)"]`)

## File to create: `msi_study/msi_study.qmd`

### Structure (one chunk per section)

**1. Setup & read data**
- Source: `code/find_many_similar.R`, `code/read_annotated_vcf.R`, `code/Generate_Koh89_Koh476_catalog_0121.R`
- `library(tidyverse)`, `library(knitr)`
- Read `Manuscript_data/Liu_et_al_476_type_spectra.tsv` → `spectra` (row.names = 1, check.names = FALSE)
- Read `Manuscript_data/Liu_et_al_final_476_type_signatures.tsv` → `sigs` (row.names = 1, check.names = FALSE)

**2. Find similar spectra**
- Call `find_many_similar()` with `sig_path`, `sig_col = params$sig_to_report`, `spectra_path`, `cosine_cutoff = 0.99`, `num_exemplars = 5`, `min_mutations = 50`, `out_pdf` based on sig name

**3. Plot exemplar spectra**
- Loop over `similar_spectra$plots` and print each plot

**4. Validate computed vs file spectra**
- For each sample ID in `similar_spectra$top_exemplars$spectrum`:
  - `tryCatch` around `read_annotated_vcf(id)`
  - If successful, add a `sample_id` column, call `GenerateKoh476CatalogfromAnnotateVcf(annot_vcf, "sample_id")`
  - Compare computed spectrum to the corresponding column in `spectra`: print total mutations (computed vs file), sum of absolute differences, cosine similarity
- `rbind` all successfully read `annot_vcf` into one combined data.table

**5. Poly-T and Poly-C tract length distributions**
- Filter combined `annot_vcf` for each of: `Del\(T`, `Ins\(T`, `Del\(C`, `Ins\(C` via `grepl()` on `Koh_476`
- For each: `count(R_outside_ins_or_del_seq)`, `mutate(repeat_count = 1 + R_outside_ins_or_del_seq)`
- Print mean, median, sd of `repeat_count`
- Plot histogram of `repeat_count`

**6. Investigate types_of_interest**
- For each `indel_type` in `params$types_of_interest`:
  - Filter `annot_vcf` where `Koh_476 == indel_type`
  - `count(ins_or_del_seq, R_outside_ins_or_del_seq, short_visual) |> arrange(desc(n))`
  - Print as kable
  - Plot histogram of `n`

## Existing functions to reuse
- `find_many_similar()` — `code/find_many_similar.R`
- `read_annotated_vcf()` — `code/read_annotated_vcf.R`
- `GenerateKoh476CatalogfromAnnotateVcf()` — `code/Generate_Koh89_Koh476_catalog_0121.R`
- `filter_vcf_by_mutation_type()` — `code/filter_vcf_by_mutation_type.R`
- Cosine similarity: use `philentropy::cosine_dist` or `lsa::cosine` (already used in codebase)

## Verification
- Render with: `quarto render msi_study/msi_study.qmd -P sig_to_report:InsDel_J -P "types_of_interest:Del2:U1:R(5,9)"`
- Check the HTML output opens and contains: exemplar plots, validation table, histograms, kable tables
