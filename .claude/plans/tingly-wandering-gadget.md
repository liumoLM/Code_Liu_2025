# Plan: Add early-stop condition to max_neg_fraction loop

## Context
In `vignette/vhelpers.R` lines 222-232, a loop increments `max_neg_fraction` from 0.005 to 0.15 searching for the highest value where `prob_ge_total_negative >= 0.5`. Currently there's no guard against the most negative residual channel becoming too large relative to the spectrum. We need to stop early if the largest negative residual exceeds 5% of the spectrum's peak.

## Change

**File:** `vignette/vhelpers.R`, lines ~222-232

Inside the `for (max_neg_fraction in seq(...))` loop, after calling `mSigBG::max_subtract_signature()`, compute the residual and break if `abs(min(residual)) > 0.05 * max(spectrum)`:

```r
for (max_neg_fraction in seq(0.005, 0.15, by = 0.005)) {
  mss_result <- mSigBG::max_subtract_signature(
    spectrum = spectrum,
    sig_to_subtract = sig_to_subtract,
    max_neg_fraction = max_neg_fraction
  )
  # Stop if most negative residual channel exceeds 5% of spectrum peak
  residual <- as.numeric(spectrum) - mss_result$n_subtract * as.numeric(sig_to_subtract)
  if (abs(min(residual)) > 0.05 * max(as.numeric(spectrum))) {
    break
  }
  if (mss_result$prob_ge_total_negative >= min_prob_ge_total_negative) {
    best_max_neg_fraction <- max_neg_fraction
    best_mss_result <- mss_result
  }
}
```

The `break` is placed before the best-result update so that the offending `max_neg_fraction` is never recorded as `best`.

## Verification
- Render the vignette and confirm the residual plots look reasonable (no excessively negative channels).
- Check the assignment log for affected signatures to confirm `max_neg_fraction` values are the same or lower than before.
