source(here::here("code_for_internal_exploration/find_many_similar.R"))

out_dir <- here::here("code_for_internal_exploration/similar_to_indesdel_2abc")
data_dir <- here::here("Manuscript_data/finalized_cap9")

for (sig_col in c("InsDel2a", "InsDel2b", "InsDel2c")) {
  for (type in c("476", "89")) {
    message("=== ", sig_col, " / ", type, "-type ===")
    res <- find_many_similar(
      sig_path       = file.path(data_dir, paste0("liu_et_al_", type, "_signatures.tsv")),
      sig_col        = sig_col,
      spectra_path   = file.path(data_dir, paste0("liu_et_al_", type, "_spectra.tsv")),
      cosine_cutoff  = 0.9,
      num_exemplars  = 300,
      do_plot        = FALSE
    )

    write.table(
      res$above_cutoff,
      file.path(out_dir, paste0("similar_to_", sig_col, "_", type, ".tsv")),
      sep = "\t", row.names = FALSE, quote = FALSE
    )

    message(sig_col, " ", type, "-type: ", nrow(res$above_cutoff), " spectra above cutoff\n")
  }
}
