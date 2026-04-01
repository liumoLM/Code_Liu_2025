# Add signature partial spectrum column to residual table

## Context
The "Per-channel residual breakdown" table in `onesig.qmd` currently shows Channel, Residual, Spectrum, and Sig_Prop. The user wants to see the per-channel counts attributed to the signature (i.e., `target_sig_partial_spectrum = n_subtract × sig_proportion`). The "residual spectrum" column already exists as the "Residual" column.

## Changes

### 1. `vignette/vhelpers.R` (~line 272-278)
Add a `Sig_Count` column to the `residual_table` data.frame:
```r
result$residual_table <- data.frame(
  Channel = ch_names,
  Spectrum = as.numeric(spectrum),
  Sig_Count = round(as.numeric(result$target_sig_partial_spectrum), 1),
  Residual = round(residual_vec, 1),
  Sig_Prop = round(as.numeric(sig_to_subtract), 5),
  stringsAsFactors = FALSE
)
```
- Reorder columns so the decomposition reads naturally: Spectrum → Sig_Count → Residual
- `Sig_Count` = per-channel mutation counts attributed to the signature

### 2. `vignette/onesig.qmd` (~line 106-113)
No changes needed — the table rendering code uses `sig_data$residual_table` directly, so the new column will appear automatically. The red-highlighting logic only touches the `Residual` column, so it remains correct.

## Verification
- Re-render the vignette and check that the "Per-channel residual breakdown" table now shows the Sig_Count column
- Verify that Spectrum ≈ Sig_Count + Residual for each row (they should sum correctly)
