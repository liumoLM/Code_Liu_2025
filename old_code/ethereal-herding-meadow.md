# Plan: Cluster CAP9 Koh89 Signatures by Cosine Similarity

## Context

67 CAP9 signature files in `Manuscript_data/Mo_CAP9_analysis/Signatures/Koh89/` each contain
columns (multinomial distributions over 89 mutation types) extracted by mSigHdp. The goal is to
cluster these ~600-700 columns by cosine similarity, find representative medoids, and filter out
small clusters — producing a consolidated set of robust signatures.

## Output

A Quarto document (`.qmd`) at the project root or in `Manuscript_data/Mo_CAP9_analysis/`.

## Implementation Steps

### 1. Load and combine signatures
- Use `rename_sigs2()` from `code/dendro2_helpers.R` to read each CAP9 file and rename columns
  (e.g., `C.PH.Bone.SoftTissue.1`, `C.H.Liver.3`)
- Assign ICAMS `catalog.row.order$ID89` as row names
- `cbind` all into a single 89-row matrix
- Track source group per column for coloring

**Key file:** `code/dendro2_helpers.R` — `rename_sigs2()` (line 34), `load_all_signatures()` (line ~170)

### 2. Compute cosine similarity matrix
- Transpose, L2-normalize rows, matrix-multiply: `cosine_sim <- norm %*% t(norm)`
- Convert to distance: `cosine_dist <- as.dist(1 - cosine_sim)`

### 3. Cluster with complete linkage
- `hc <- hclust(cosine_dist, method = "complete")`
- `clusters <- cutree(hc, h = 1 - min_similarity)` (default `min_similarity = 0.95`, cut at 0.05)
- Complete linkage guarantees all pairs within a cluster have cosine sim ≥ threshold

### 4. Find medoids
- For each cluster: medoid = member with highest average cosine similarity to other members
- `diag(sub_sim) <- NA; names(which.max(rowMeans(sub_sim, na.rm = TRUE)))`

### 5. Filter and summarize
- Drop clusters with fewer than `min_cluster_size` (default 3) members
- Summary table per retained cluster: ID, size, medoid name, min/mean/max within-cluster cosine sim, member list

### 6. Visualizations

1. **Dendrogram** — colored by source group, dashed line at cut height (0.05). Labels suppressed
   if too many leaves; color-coded instead.
2. **Cosine similarity heatmap** — ordered by dendrogram, cluster boundaries marked
3. **Silhouette plot** — `cluster::silhouette()` with cosine distance for retained clusters
4. **Medoid signature plots** — one `mSigPlot::plot_89()` per retained cluster, titled with
   cluster ID, size, and min within-cluster cosine sim. For small clusters, overlay all members
   as semi-transparent lines.
5. **Summary table** — rendered with `knitr::kable()` or `DT::datatable()`

### 7. Quality assessment (printed inline)
- Total signatures loaded, number of clusters, retained vs dropped
- Global min within-cluster cosine similarity (verify ≥ threshold)
- Silhouette width summary

### 8. Save outputs
- `*_summary.tsv` — cluster summary table
- `*_medoid_signatures.tsv` — 89-row × N-cluster medoid matrix
- `*_all_cluster_assignments.tsv` — signature → cluster mapping with retained flag

## Key files to modify/create
- **Create:** `Manuscript_data/Mo_CAP9_analysis/cluster_cap9_koh89.qmd`
- **Reuse:** `code/dendro2_helpers.R` (rename_sigs2, load pattern)
- **Reuse patterns from:** `code/cluster_catalogs.R` (dendrogram viz), `code/all_pairwise_cosine.R` (heatmap)

## Parameters (set in first code chunk)
- `min_similarity <- 0.95`
- `min_cluster_size <- 3`
- `sig_dir <- here::here("Manuscript_data", "Mo_CAP9_analysis", "Signatures", "Koh89")`

## Dependencies
`here`, `lsa` or matrix algebra, `ggplot2`, `ggdendro`, `mSigPlot`, `cluster`, `gridExtra`, `DT`

## Verification
1. Render the .qmd: `quarto render cluster_cap9_koh89.qmd`
2. Check that the global min within-cluster cosine sim is ≥ 0.95
3. Inspect medoid plots for biological plausibility
4. Review silhouette plot for any negative silhouette widths (mis-clustered signatures)
