# Plan: Add Del3/Del4 tables to study_sigF.qmd

## Context

The user wants to extend `sigFstudy/study_sigF.qmd` to show tables for additional deletion types beyond the existing Del2:M1 and Del2:U2:R2. The new types are: **Del3:U3:R2, Del4:U4:R2, Del3:M2, Del4:M3, Del4:M2**. All five types exist in the data. The tumor CSV files need regeneration (existing files have `prev_` prefix but code expects no prefix).

## Key findings

1. **CSV files need regeneration**: Files exist as `prev_DRUP01030028T_annotvcf.csv` and `prev_CPCT02100089T_annotvcf.csv` but code expects names without `prev_`. The code already handles regeneration from raw VCFs when CSVs don't exist.

2. **`del2_counts` helper needs generalization**: The `nchar(ins_or_del_seq) >= 2` filter is broad enough, but the regex `(?<=<).{2}(?=>)` only captures exactly 2 characters. For Del3/Del4 types it needs to capture variable-length sequences.

3. **`revc` function**: Comes from `ICAMS::revc` (exported, vectorized). Currently called as bare `revc()` — works because `library(ICAMS)` is loaded at the top.

4. **`short_visual` format examples**:
   - U-types: `<AGG>[AGG]` (Del3:U3:R2), `<CTTT>[CTTT]` (Del4:U4:R2)
   - M-types: `<{AG}G>{AG}` (Del3:M2), `<{TTT}C>{TTT}` (Del4:M3), `<{AC}TC>{AC}` (Del4:M2)

## Changes to `sigFstudy/study_sigF.qmd`

### Step 1: Generalize `del2_counts` → `del_counts`

Rename to `del_counts` and fix the dinuc regex to capture variable-length sequences:

```r
del_counts <- function(annotvcf) {
  annotvcf %>%
    dplyr::filter(nchar(ins_or_del_seq) >= 2) %>%
    dplyr::filter(ins_or_del == "d") %>%
    dplyr::count(short_visual, Koh_476) %>%
    dplyr::arrange(desc(n)) %>%
    dplyr::mutate(vis2 = gsub("\\}", "", gsub("\\{", "", short_visual))) %>%
    dplyr::mutate(del_seq = str_extract(vis2, "(?<=<)[^>]+(?=>)")) %>%
    dplyr::mutate(hasmh = grepl("\\{", short_visual))
}
```

- Regex: `(?<=<).{2}(?=>)` → `(?<=<)[^>]+(?=>)` (variable-length)
- Rename `dinuc` → `del_seq` (since it's no longer always a dinucleotide)

### Step 2: Update `normit` to use `del_seq` instead of `dinuc`

```r
normit <- function(xx) {
  xx %>%
    mutate(post = gsub(".*>", "", vis2)) %>%
    mutate(post = gsub("\\[", "", post)) %>%
    mutate(clean = gsub("]", "", gsub("[><\\[]", "", vis2))) %>%
    mutate(
      norm = ifelse(substr(clean, 1, 1) %in% c("A", "G"), revc(clean), clean)
    ) %>%
    mutate(prop_n = n / sum(xx$n)) %>%
    select(Koh_476, norm, n, prop_n, Koh_476, short_visual)
}
```

No functional changes needed here — `normit` uses `clean` (derived from `vis2`), not `dinuc`. The strand normalization logic works for any length.

### Step 3: Update callers of `del2_counts` → `del_counts`

Lines 106-124: Replace `del2_counts(...)` with `del_counts(...)`.

### Step 4: Add new deletion types to the table loop

Extend the `for` loop (line 141) to include the five new types:

```r
for (deltype in c(
  "Del2:M1",
  "Del2:U2:R2",
  "Del3:U3:R2",
  "Del4:U4:R2",
  "Del3:M2",
  "Del4:M3",
  "Del4:M2"
)) {
```

### Step 5: Regenerate the tumor annotated CSVs

Delete the `prev_`-prefixed CSVs so the existing code regenerates them fresh:
- `rm prev_DRUP01030028T_annotvcf.csv prev_CPCT02100089T_annotvcf.csv`

This ensures the CSVs match what the code expects and are generated with the current version of ICAMS.

## Files to modify

- `sigFstudy/study_sigF.qmd` — all changes above

## Files to delete

- `sigFstudy/prev_DRUP01030028T_annotvcf.csv`
- `sigFstudy/prev_CPCT02100089T_annotvcf.csv`

## Verification

1. Render the qmd: `cd sigFstudy && quarto render study_sigF.qmd`
2. Open the HTML and verify:
   - CSV regeneration succeeds (no errors in VCF annotation chunks)
   - All 7 deletion type tables render with data from all 5 sources
   - The `norm` column shows strand-normalized sequences of appropriate length (2, 3, or 4 bases)
