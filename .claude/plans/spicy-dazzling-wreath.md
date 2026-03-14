# Plan: Show top 3 best-matching spectra per medoid

## Context
The current code (lines 596–611 of `cluster_cap9_koh89.qmd`) plots only the single best-matching spectrum for each medoid. The user wants to see the top 3 best-matching spectra instead.

## File to modify
- `Manuscript_data/Mo_CAP9_analysis/cluster_cap9_koh89.qmd` — the `medoid-plots` chunk (lines 596–611)

## Changes
Replace the current single-best-match block with a loop over the top 3:

1. Sort `cos_sims` in decreasing order and take the first 3.
2. Loop over these 3, printing the spectrum name + cosine similarity and plotting each with `plot_89()`.
3. Label them "Best-matching spectrum 1/2/3" in the plot title.

## Verification
- Render with `Rscript Manuscript_data/Mo_CAP9_analysis/render_cluster.R`
- Check that each medoid section now shows 3 best-matching spectrum plots
