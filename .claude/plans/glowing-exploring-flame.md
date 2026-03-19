# Plan: Set Jin et al minimum cosine similarity to 0.9

## Context
The vignette currently shows all Jin et al. signature matches regardless of cosine similarity (`Jin_min_cosine <- 0.`). This clutters the output with poor matches. The user wants only matches with cosine >= 0.9, consistent with the COSMIC and Koh thresholds already set to 0.9.

## Change

### File: `vignette/vignette.qmd` (line 69)

Change:
```r
Jin_min_cosine <- 0.
```
to:
```r
Jin_min_cosine <- 0.9
```

No other files need changes — this threshold already controls:
- Which Jin matches are stored in `jin_matches` (lines 228-241)
- Which Jin plots are generated per signature (`vhelpers.R` lines 758-796)
- Which Jin matches appear in the overview table (lines 479-533)
- Which Jin plots are shown in `onesig.qmd` (lines 407-434)

## Verification
- Render the vignette and confirm Jin matches only appear where cosine >= 0.9
- Check the overview table's "Closest Sig from Jin et al." column has NA for signatures without a 0.9+ match
