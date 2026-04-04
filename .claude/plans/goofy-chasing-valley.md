# Add nonzero/total ratio to hamburger plot sample count annotations

## Context
The hamburger plots already show two numbers below each cancer type column: the number of samples with non-zero signature activity (top) and the total number of samples (bottom). The user wants to also display the ratio (nonzero / total) so the proportion is immediately visible without mental math.

## Changes

### 1. `code_for_internal_exploration/hamburger_plots/hamburger_interactive.qmd` (lines 299–313)

Currently the plotly annotation text is:
```r
text = paste0(count_df$nonzero_samples[i], "<br>", count_df$total_samples[i])
```

Change to include the ratio as a third line:
```r
text = paste0(
  count_df$nonzero_samples[i], "<br>",
  count_df$total_samples[i], "<br>",
  round(count_df$nonzero_samples[i] / count_df$total_samples[i], 2)
)
```

Also increase bottom margin slightly (b = 70 → ~85) to accommodate the extra line.

### 2. `code_for_internal_exploration/plot_signature_assignments.R` (lines 306–317)

Currently the ggplot2 `geom_text` label is:
```r
label = paste0(nonzero_samples, "\n", total_samples)
```

Change to include the ratio as a third line:
```r
label = paste0(nonzero_samples, "\n", total_samples, "\n",
               round(nonzero_samples / total_samples, 2))
```

Also adjust `vjust` or `plot.margin` bottom if needed to prevent clipping of the extra line.

## Verification
- Run the Shiny app: `quarto serve hamburger_interactive.qmd` and confirm three lines appear below each cancer type column (nonzero, total, ratio).
- Run `Rscript plot_signature_assignments.R` and open the output PDF to confirm the ratio line appears below each column.
