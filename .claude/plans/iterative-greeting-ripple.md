# Plan: Add logging for `assignment` and `sigid` in vhelpers.R

## Context
The user wants to debug the residual spectrum calculation near line 250-253 of `vignette/vhelpers.R` by logging the values of `assignment` (including row names) and `sigid` to a log file.

## File to modify
- `vignette/vhelpers.R` (line ~253, after `result$residual_spectrum` is computed)

## Change
Insert logging code after line 253 (`result$residual_spectrum[result$residual_spectrum < 0] <- 0`) that:

1. Opens/appends to `vignette/assign_for_residual.log`
2. Writes `sigid` value
3. Writes `assignment` matrix with row names
4. Closes the connection

```r
        # Log assignment and sigid for debugging
        log_con <- file(here::here("vignette", "assign_for_residual.log"), open = "a")
        writeLines(paste0("=== sigid: ", sigid, " ==="), log_con)
        writeLines("assignment (with row names):", log_con)
        capture.output(print(assignment), file = log_con, append = TRUE)
        writeLines("", log_con)
        close(log_con)
```

This goes between line 253 and the `if (FALSE) {` block on line 255.

## Verification
- Render the vignette that calls this function and check that `vignette/assign_for_residual.log` is created with the expected content.
