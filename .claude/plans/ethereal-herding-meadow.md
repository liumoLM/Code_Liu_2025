# Plan: Add medoid pairwise cosine similarity heatmap to cluster_cap9_koh89.qmd

## Context
The file `code/all_pairwise_cosine.R` generates a lower-triangle ggplot2 heatmap of pairwise cosine similarities between signatures (e.g., `plot_output/heatmap_89_type_signatures.pdf`). The user wants a similar heatmap added to `cluster_cap9_koh89.qmd`, but using the **medoid signatures** from retained clusters instead of the full signature set.

## Changes

### Modify: `Manuscript_data/Mo_CAP9_analysis/cluster_cap9_koh89.qmd`

Add a new chunk after the dendrogram section (before "Summary statistics"), roughly following the approach in `code/all_pairwise_cosine.R` (lines 44–131):

1. **Compute medoid pairwise cosine similarities** using `cosine_sim` (already computed) subsetted to `names(medoid_names)`.
2. **Build lower-triangle heatmap** with ggplot2:
   - `geom_tile()` colored white→red
   - `geom_text()` labels for cells ≥ 0.9
   - Vertical line segments from high-similarity cells to x-axis (as in the original)
   - Axis labels = medoid names (stripped of "C." prefix to match dendrogram)
3. Use `reshape2::melt()` for the matrix→long conversion (already a dependency via the existing code pattern).

**Available variables at insertion point:**
- `cosine_sim` — full pairwise cosine similarity matrix (all signatures)
- `medoid_names` — named vector: names = medoid sig names, values = cluster IDs
- `combined` — full signature matrix

**Libraries needed:** `reshape2` (add to setup chunk).

### No changes to `code/dendro2_helpers.R` or `code/all_pairwise_cosine.R`

The heatmap code will be inline in the qmd chunk, adapted from `all_pairwise_cosine.R` lines 60–131.

## Verification
`quarto render cluster_cap9_koh89.qmd` and visually check the heatmap shows medoid-vs-medoid cosine similarities.
