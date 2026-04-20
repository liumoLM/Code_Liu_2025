# Plan: `fig1.R` — 3-panel stacked signature plot (portrait letter)

## Context

The user wants a first-draft Figure 1 script for the paper. It stacks three rows — one each for `plot_ID83`, `plot_ID89`, `plot_ID476` — on a portrait letter PDF, showing a single biologically linked signature (`C_ID23` ↔ `InsDel23`) under all three classifications. The plots should be small (roughly the top 1/3–1/2 of the page) so the rest of the page is free for annotation in a PDF editor (Illustrator/Inkscape). PDF text must remain editable.

## Decisions (from user)

- Signatures: `C_ID23` in 83, `InsDel23` in 89 and 476 (biologically matched per `connection_table.tsv`)
- Page: **portrait** letter (8.5 × 11 in), plots span printable width **6.5 in** (1 in side margins)
- Arrows: `num_peak_labels = 4` on all three plots
- Plot block vertical: ~4.5 in tall total (middle of the user's 1/3–1/2 range); stays near the top of the page; rest blank

## Inputs

Read each TSV with `read.delim(..., row.names = 1, check.names = FALSE)`. Confirmed these paths and column names exist:

- `Manuscript_data/finalized_cap9/liu_et_al_83_signatures.tsv` — col `C_ID23`
- `Manuscript_data/finalized_cap9/liu_et_al_89_signatures.tsv` — col `InsDel23`
- `Manuscript_data/finalized_cap9/liu_et_al_476_signatures.tsv` — col `InsDel23`

Pass each to the respective plot function as a **single-column data frame** (`sigs[, col, drop = FALSE]`) — this is the idiom already used in `code_for_internal_exploration/plot_all_signatures.R:33-46`.

## Output

`/home/steve/github/Code_Liu_2025/some_figures/fig1.pdf` (written by `fig1.R` in the same directory).

## Approach

New file: `/home/steve/github/Code_Liu_2025/some_figures/fig1.R`

Structure:

1. **Header + params block** — top-of-file variables so the user can tweak without digging:
   - `sig_col_83 <- "C_ID23"`, `sig_col_89 <- "InsDel23"`, `sig_col_476 <- "InsDel23"`
   - `num_peaks <- 4`
   - `base_size_83 <- 11`, `base_size_89 <- 11`, `base_size_476 <- 11` — per-plot font scaling, each passed to its matching `plot_ID*(..., base_size = ...)` call
   - `page_w <- 8.5`, `page_h <- 11` (portrait letter)
   - `plot_w <- 6.5`, `plot_h <- 4.5` (plot-block dimensions; edit to resize)
   - `margin_left <- 1`, `margin_top <- 1`
   - `out_file <- here::here("some_figures/fig1.pdf")`

2. **Load packages**: `ggplot2`, `patchwork`, `grid`, `here`, `mSigPlot`. Use `mSigPlot::` prefix rather than `library(mSigPlot)`.

3. **Read inputs** with `here::here("Manuscript_data/finalized_cap9/...")` and `read.delim(..., row.names = 1, check.names = FALSE)`.

4. **Build the three ggplots** (each returns a ggplot); pass the matching `base_size_*` so each panel's font size is independently tunable:
   ```r
   p83  <- mSigPlot::plot_ID83 (sigs_83 [, sig_col_83,  drop = FALSE], num_peak_labels = num_peaks, base_size = base_size_83)
   p89  <- mSigPlot::plot_ID89 (sigs_89 [, sig_col_89,  drop = FALSE], num_peak_labels = num_peaks, base_size = base_size_89)
   p476 <- mSigPlot::plot_ID476(sigs_476[, sig_col_476, drop = FALSE], num_peak_labels = num_peaks, base_size = base_size_476)
   ```

5. **Stack vertically** with `patchwork`: `combined <- p83 / p89 / p476`.

6. **Position inside the page** using `cairo_pdf()` + a `grid::viewport` so the plot block sits in the top-left printable area (1 in from left, 1 in from top, 6.5 in wide, 4.5 in tall):
   ```r
   cairo_pdf(out_file, width = page_w, height = page_h)
   grid::grid.newpage()
   vp <- grid::viewport(
     x      = grid::unit(margin_left, "in"),
     y      = grid::unit(page_h - margin_top, "in"),
     width  = grid::unit(plot_w, "in"),
     height = grid::unit(plot_h, "in"),
     just   = c("left", "top"))
   grid::pushViewport(vp)
   print(combined, newpage = FALSE)
   grid::popViewport()
   dev.off()
   ```

   `cairo_pdf()` is what keeps text editable in Illustrator/Inkscape (same idiom already used in `vignette/old/test_473_to_83.R:33` and `code/all_pairwise_89_476.R:257-262`). No raster fallbacks; no `useDingbats` flag needed.

## Things NOT to do

- Do not use `grDevices::pdf()` — fonts/text may get outlined depending on settings; `cairo_pdf` is the safer default for editable text.
- Do not hard-code dimensions inside plot function calls; keep all sizing in the params block so "we'll adjust later" stays cheap.
- Do not call `library(dplyr)` or similar; no `dplyr::filter` needed in this script.

## Critical files

- **New**: `/home/steve/github/Code_Liu_2025/some_figures/fig1.R`
- **Reads**: the three TSVs in `Manuscript_data/finalized_cap9/`
- **Writes**: `/home/steve/github/Code_Liu_2025/some_figures/fig1.pdf`
- **Reference idioms**: `code_for_internal_exploration/plot_all_signatures.R:33-46` (TSV → plot), `vignette/old/test_473_to_83.R:33` (`cairo_pdf` usage)

## Verification

1. Run `Rscript some_figures/fig1.R` from the repo root; expect no warnings other than ggrepel overlap messages.
2. `xdg-open some_figures/fig1.pdf` and check:
   - Three panels stacked in the top portion of a portrait letter page.
   - Each panel shows peak-label arrows (4 per plot).
   - Text in the PDF is selectable (not outlined) — quick test: open in a PDF viewer and try to select any axis label.
3. Sanity-check the signatures: the three panels should look like the same biological signature at three resolutions (they are related per `connection_table.tsv`).
4. If the plot block feels too tall/short, edit `plot_h` (and optionally `plot_w`) in the params block and re-run. Tune each panel's text separately via `base_size_83` / `base_size_89` / `base_size_476`.
