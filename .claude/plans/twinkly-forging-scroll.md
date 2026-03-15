# Plan: Show best-matching converted-476 signature per cluster

## Context
When `mutation_type == 89`, the 476-type medoid signatures converted to 89-type are already loaded as `cv476_sigs` in the `medoid-ref-load` chunk. The user wants each cluster's medoid plot section to also show the most similar signature from this converted set.

## Change
**File:** `Manuscript_data/Mo_CAP9_analysis/cluster_cap9.qmd`, lines ~632-653 (in the `medoid-plots` chunk)

After the best-matching spectra loop (line 652), add a block guarded by `if (mutation_type == 89)` that:
1. Computes cosine similarity between the medoid and each column of `cv476_sigs`
2. Finds the best match
3. Prints a markdown header and the plot using `plot_89()`

## Verification
Render with: `Rscript Manuscript_data/Mo_CAP9_analysis/render_cluster_cap9.R --min-similarity 0.95`
Check that the 89-type output shows the converted-476 best match after each cluster's spectra section.
