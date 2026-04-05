# Add microhomology length column to indel_deep_dive output

## Context
The `short_visual` field uses curly braces `{...}` to denote the microhomology sequence in MH-type deletions (e.g., `<{ACTCGAA}G>{ACTCGAA}`). We want a new column `mh_length` = number of characters between `{` and `}` (NA if no curly braces).

## File to modify
- `code_for_internal_exploration/indel_deep_dive/indel_deep_dive.R`

## Change
After computing `R_intuitive`, add:
```r
visual_table$mh_length <- ifelse(
  grepl("\\{", visual_table$short_visual),
  nchar(sub(".*\\{([^}]+)\\}.*", "\\1", visual_table$short_visual)),
  NA_integer_
)
```

Add `mh_length` to the `dplyr::select()` and empty tibble definition.

## Verification
- Run `Rscript deep_dive_H.R` and check that MH rows have integer `mh_length`, repeat rows have NA.
