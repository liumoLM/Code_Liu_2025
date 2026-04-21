# deep_dive_for_fig1.R — one-off
#
# For every 476-type spectrum with cosine similarity > 0.90 to signature
# InsDel23, run indel_deep_dive() keyed on the top peaks of InsDel23, then
# sum the per-tumor counts into a single aggregated table and save it to
# some_figures/fig1_deep_dive.csv.

library(here)
library(dplyr)

source(here::here("code_for_internal_exploration/indel_deep_dive.R"))

sig_path     <- here::here(
  "Manuscript_data/finalized_cap9/liu_et_al_476_signatures.tsv")
spectra_path <- here::here(
  "Manuscript_data/finalized_cap9/liu_et_al_476_spectra.tsv")
sig_col        <- "InsDel23"
cosine_cutoff  <- 0.90

# ---- Samples: cosine > cosine_cutoff to InsDel23 --------------------------
# Inlined (avoids find_many_similar(), which calls a stale plot_476()).

sigs    <- read.delim(sig_path,     row.names = 1, check.names = FALSE)
spectra <- read.delim(spectra_path, row.names = 1, check.names = FALSE)
stopifnot(nrow(sigs) == nrow(spectra))

s <- sigs[, sig_col]
s <- s / sum(s)

cossim <- function(a, b) sum(a * b) / sqrt(sum(a * a) * sum(b * b))
cos_vals <- vapply(seq_len(ncol(spectra)), function(j) {
  v <- as.numeric(spectra[, j])
  if (sum(v) == 0) return(NA_real_)
  cossim(s, v / sum(v))
}, numeric(1))
names(cos_vals) <- colnames(spectra)

samples <- names(sort(cos_vals[cos_vals > cosine_cutoff], decreasing = TRUE))
cat(sprintf("Spectra with cosine > %.2f to %s: %d\n",
            cosine_cutoff, sig_col, length(samples)))
cat("  ", paste0(samples, " (", sprintf("%.4f", cos_vals[samples]), ")",
                 collapse = ", "), "\n", sep = "")

# ---- Patterns: InsDel23 peaks contributing >= 1% --------------------------

sigs    <- read.delim(sig_path, row.names = 1, check.names = FALSE)
sig_vec <- setNames(sigs[, sig_col], rownames(sigs))

# Escape regex metacharacters so signature rownames match literally in
# indel_deep_dive()'s grepl() call against the Koh_476 column.
escape_re <- function(s) {
  meta <- c("\\", "[", "]", "(", ")", ".", "+", "*", "?", "^", "$", "|",
            "{", "}")
  for (ch in meta) s <- gsub(ch, paste0("\\", ch), s, fixed = TRUE)
  s
}

peak_names <- names(sig_vec)[sig_vec >= 0.01]
peak_names <- peak_names[order(-sig_vec[peak_names])]
cat(sprintf("InsDel23 peaks >= 1%% each: %d (total mass = %.3f)\n",
            length(peak_names), sum(sig_vec[peak_names])))

patterns_to_match <- setNames(
  peak_names,
  paste0("^", escape_re(peak_names), "$"))

# ---- Run indel_deep_dive -------------------------------------------------

sample_info <- read.delim(
  here::here("Manuscript_data/sample_info.tsv"), check.names = FALSE)

pdf_path <- tempfile("deep_dive_for_fig1_", fileext = ".pdf")
visual_table <- indel_deep_dive(
  samples_to_fetch  = samples,
  patterns_to_match = patterns_to_match,
  sample_info       = sample_info,
  cap9              = FALSE,
  pdf_path          = pdf_path)

# ---- Sum matching rows across tumors -------------------------------------

summed <- visual_table |>
  dplyr::group_by(rc_visual, visual, pattern, short_visual, ins_or_del_seq,
                  R, R_intuitive, L, mh_length, prefix, suffix) |>
  dplyr::summarise(count = sum(count), .groups = "drop") |>
  dplyr::select(rc_visual, visual, count, pattern, short_visual,
                ins_or_del_seq, R, R_intuitive, L, mh_length,
                prefix, suffix) |>
  dplyr::arrange(dplyr::desc(count))

prelim_csv <- here::here("some_figures", "deep_dive_for_fig1_prelim.csv")
write.csv(summed, prelim_csv, row.names = FALSE)
cat(sprintf("Wrote %s (%d rows)\n", prelim_csv, nrow(summed)))

# Second output: collapse the prelim on rc_visual alone, summing counts.
summed_rc <- summed |>
  dplyr::group_by(rc_visual, pattern) |>
  dplyr::summarise(count = sum(count), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(count))

out_csv <- here::here("some_figures", "deep_dive_for_fig1.csv")
write.csv(summed_rc, out_csv, row.names = FALSE)
cat(sprintf("Wrote %s (%d rows)\n", out_csv, nrow(summed_rc)))
