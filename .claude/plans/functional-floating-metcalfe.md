# Plan: Dendrogram improvements (height, label filtering, cluster dropdown)

## Context

The dendrograms in `cluster_cap9_koh89.qmd` have too many leaf labels (~621 signatures), making them unreadable. The user wants: (1) 1.5x more vertical space, (2) fewer leaf labels, and (3) a dropdown to jump to specific cluster/medoid positions.

## Files to modify

- `code/dendro2_helpers.R` — `build_plotly_dendrogram()` function (lines 577-711)
- `Manuscript_data/Mo_CAP9_analysis/cluster_cap9_koh89.qmd` — three dendrogram blocks

## Changes

### 1. Add `label_filter` parameter to `build_plotly_dendrogram()`

New parameter: `label_filter = NULL` (character vector of leaf names to annotate; NULL = all).

In the annotations loop (line 644), filter to only create text annotations for leaves in `label_filter`. Markers (dots) and hover remain for ALL leaves.

### 2. Add `cluster_dropdown` parameter to `build_plotly_dendrogram()`

New parameter: `cluster_dropdown = FALSE`.

When TRUE and `medoids`/`clusters` are provided:
- Build plotly `updatemenus` dropdown with one button per cluster: `"<clusterID> <medoidName>"`
- Each button calls `relayout` to set `xaxis.range` centered on that cluster's member span + padding
- Include a "Show All" button at the top
- Sort buttons by cluster ID
- Position dropdown at top-left (x=0, y=1.15); increase top margin to 60px when dropdown is active

### 3. Increase heights in QMD

| Block | Current | New |
|-------|---------|-----|
| dendrogram (line 214) | no explicit height | `\|> layout(height = 1350)` |
| medoid-ref-dendrogram-med (line 310) | `layout(height = 900)` | `layout(height = 1350)` |
| medoid-ref-dendrogram-all (line 338) | `layout(height = 900)` | `layout(height = 1350)` |

### 4. Pass new parameters in QMD

| Block | `label_filter` | `cluster_dropdown` |
|-------|---------------|-------------------|
| dendrogram (~621 leaves) | `names(medoid_names)` | `TRUE` |
| medoid-ref-dendrogram-med (few leaves) | omit (show all) | `TRUE` |
| medoid-ref-dendrogram-all (~621+ leaves) | `names(medoid_names)` | `TRUE` |

## Backward compatibility

Both new parameters default to preserving old behavior (NULL / FALSE), so `interactive_dendro2.qmd` is unaffected.

## Verification

Render with: `Rscript Manuscript_data/Mo_CAP9_analysis/render_cluster.R --min-similarity 0.95 --mutation-type 89`

Check:
- Dendrograms are taller
- Only medoid labels shown (with cluster ID prefix) on the all-signatures dendrograms
- Dropdown navigates to correct cluster positions
- Medoid-ref dendrogram still shows all labels
