# Plan: Show reference signature labels in medoid-ref-dendrogram-all

## Context

In the "Medoid + reference dendrogram" section of `cluster_cap9_koh89.qmd`, the `medoid-ref-dendrogram-all` block currently sets `label_filter = names(medoid_names)`, so only medoid extraction signatures get text labels. Reference signatures (Liu, Koh) are present as dots with hover but have no visible text labels. The user wants reference signature names shown as well.

## File to modify

- `Manuscript_data/Mo_CAP9_analysis/cluster_cap9_koh89.qmd` — `medoid-ref-dendrogram-all` block (line ~334)

## Change

In the `medoid-ref-dendrogram-all` block, expand `label_filter` to include reference signature names alongside medoid names:

```r
label_filter = c(names(medoid_names), colnames(liu_sigs),
                 if (mutation_type == 89) colnames(koh_sigs))
```

This uses variables already in scope from the `medoid-ref-load` block (line 237): `liu_sigs`, `koh_sigs`, `mutation_type`.

No changes needed to `build_plotly_dendrogram()` — the `label_filter` parameter already accepts any character vector of leaf names.

## Verification

Render with: `Rscript Manuscript_data/Mo_CAP9_analysis/render_cluster.R --min-similarity 0.95 --mutation-type 89`

Check that the `medoid-ref-dendrogram-all` dendrogram shows text labels for Liu sigs, Koh sigs, and medoids, but NOT for the ~600 non-medoid extraction signatures.
