## MSIStudy

Specification for the study of MSI tumors, 

- Reports on the length distributions of poly-T tracts with T insertions or deletions
- Reports on the length distributions of poly-T tracts with T insertions or deletions
- Reports on sequences and lengths involved in deletions or insertions of > 1 base
- This should be in a .qmd document with the user arguments supplied in whatever way is common. This arguments are sig_to_report and types_of_interest

## Inputs

- Manuscript_data/Liu_et_al_476_type_spectra.tsv
- Manuscript_data/Liu_et_al_final_476_type_signatures.tsv
- The id of indel signature; call this sig_to_report
    (usually one of "InsDel2b", "InsDel2c", "InsDel7", "InsDel_J", "InsDel_Kbeta", "InsDel_O", "InsDel_P")
- A list of 'types_of_interest' which are a subset of rownames Manuscript_data/Liu_et_al_final_476_type_signatures.tsv
  Defaults to "Del2:U1:R(5,9)"

-- source code/find_many_similar.R
-- source code/read_annoated_vcf.R
-- source code/Generate_Koh89_Koh476_catalog_0121.R

## Ouput

The code will be quarto markdown, and the ouput will be the result of rendering the .qmd.

## Algorithm

Write R code. Use tidyverse.

### Read file 
Read Manuscript_data/Liu_et_al_476_type_spectra.tsv into variable spectra; the first column row names; columns are indexed by sample id, rows are indexed by mutation type

Read Manuscript_data/Liu_et_al_final_476_type_signatures.tsv into variable sigs; the first column is row names; columns are indexed by signature_id, rows are indexed by mutation type (same order as the rows of spectra)

### Find spectra similar to the sig_to_report
  call this (R) code, 
  
   ```{r}
similar_spectra <- find_many_similar(
  sig_path = here::here(
    "Manuscript_data/Liu_et_al_final_476_type_signatures.tsv"
  ),
  signature_id,
  here::here("Manuscript_data/Liu_et_al_476_type_spectra.tsv"),
  cosine_cutoff = 0.99
  num_exemplars = 5,
  out_pdf = "exemplars_for_IndDelJ.pdf",
  min_mutations = 50
)
```

### Plot the spectra you found

Put all the plots in similar_spectra$plot in the rendered quarto output.

### Check newly computed spectra to spectra read from from the file
For id (column name) I think in similar_spectra do:

    annot_vcf <- read_annotated_vcf(id)
    use a do-try; if the it doesn't work go on to the next

    If you get the vcf, call GenerateKoh476CatalogfromAnnotateVcf()
    on it to get a computed_spectrum.

    print out the total number of mutations in the computed_spectrum, land
    the number of mutations in the corresonding
    real spectrum from similar_spectra.
    print out the sum of the absolute differences between computed
    spectrum and the real specgrum
    print out the cosine similarity between the two


rbind all the annot_vcf can call it annot_vcf

### Report on the length distributions of poly-T tracts with T insertions or deletions and similar

You can use this code to select the row for deletions of single Ts:

polyt <- annot_vcf |> dplyr::filter(grepl("Del\\(T", Koh_476)) |>
    dplyr::count(R_outside_ins_or_del_seq)

Write analogous code fir Ins\\(T,  Del\\(C, Ins\\(C

in the output polyt (or whaever, mutate(repeat_count = 1 + R_outside_ins_or_del_seq))

write mean, meadian, sd, and plot a histogram of repeat_count

### Investigate the types_of_interest

for each indel_type in types_of_interest do:

```{r}
annot_vcf |> filter(Koh_476 == indel_type) |>
  dplyr::count(ins_or_del_seq, R_outside_ins_or_del_seq, short_visual) |>
  dplyr::arrange(desc(n)) -> summary_of_type
```
Print summary_of_type as a kagle.

Plot a histogram of n.
  
