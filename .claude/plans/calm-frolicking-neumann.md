# Plan: Add `max_neg_fraction` to results and display

## Context
The `compute_sig_data()` function in `vhelpers.R` calls `max_subtract_signature` with a `max_neg_fraction` argument (initially 0.02, possibly bumped to 0.05). This value is not currently recorded in `result`, not written to `assign_for_residual.md`, and not shown in the vignette fold. The user wants it tracked and displayed throughout.

## Changes

### 1. `vignette/vhelpers.R` — Record `max_neg_fraction` in result

**Line ~252** (after `result$prob_ge_total_negative`): Add:
```r
result$max_neg_fraction <- max_neg_fraction
```

**Line ~359** (else branch for InsDel15/16): Add:
```r
result$max_neg_fraction <- NA
```

### 2. `vignette/vhelpers.R` — Write `max_neg_fraction` to `assign_for_residual.md`

**Line ~318** (in the markdown logging block, after the `## type89_sig_id` header and before `N subtract`): Add:
```r
writeLines(sprintf("- Max neg fraction: %.2f", max_neg_fraction), log_con)
```

### 3. `vignette/onesig.qmd` — Display with newlines and descriptions

**Lines ~82-95**: Replace the single-line `|`-separated `sprintf` with a multi-line version using `<br>` (since this is inside an HTML `<details>` block), with a short description for each value:

```r
cat(sprintf(
  paste0(
    "\nMax neg fraction (maximum allowed negative residual as fraction of total): %.2f<br>",
    "N subtract (mutation count assigned to signature): %.1f<br>",
    "N residual (mutation count in residual): %.1f<br>",
    "Total negative (sum of negative residual channels): %.1f<br>",
    "N negative channels (number of channels with negative residual): %d<br>",
    "P(\u2265 total negative) (probability of seeing this much negativity by chance): %.3f\n\n"
  ),
  sig_data$max_neg_fraction,
  sig_data$n_subtract,
  sig_data$n_residual,
  sig_data$total_negative,
  sig_data$n_negative_channels,
  sig_data$prob_ge_total_negative
))
```

### 4. `vignette/vignette.qmd` — Exclude from `sig_table2`

Already excluded: line 498 has `-n_subtract, -n_residual, -total_negative, -n_negative_channels, -prob_ge_total_negative, -residual_table`. Add `-max_neg_fraction` to this list.

## Files to modify
- `vignette/vhelpers.R` (3 insertions)
- `vignette/onesig.qmd` (replace sprintf block)
- `vignette/vignette.qmd` (add `-max_neg_fraction` to select exclusion)

## Verification
- Re-render the vignette: `cd vignette && Rscript render_cluster.R` (or appropriate render command)
- Check `assign_for_residual.md` for the new `Max neg fraction` line under each signature
- Check the HTML output's fold sections for the new multi-line format with descriptions
