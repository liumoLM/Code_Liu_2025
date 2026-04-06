source(here::here("code_for_internal_exploration/find_many_similar.R"))

out_dir <- here::here("code_for_internal_exploration/similar_to_indesdel_N_beta")
data_dir <- here::here("Manuscript_data/finalized_cap9")

# 476-type analysis
res_476 <- find_many_similar(
  sig_path     = file.path(data_dir, "liu_et_al_476_signatures.tsv"),
  sig_col      = "InsDel_N_beta",
  spectra_path = file.path(data_dir, "liu_et_al_476_spectra.tsv"),
  cosine_cutoff  = 0.9,
  num_exemplars  = 300,
  do_plot        = FALSE
)

write.table(
  res_476$above_cutoff,
  file.path(out_dir, "similar_to_InsDel_N_beta_476.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# 89-type analysis
res_89 <- find_many_similar(
  sig_path     = file.path(data_dir, "liu_et_al_89_signatures.tsv"),
  sig_col      = "InsDel_N_beta",
  spectra_path = file.path(data_dir, "liu_et_al_89_spectra.tsv"),
  cosine_cutoff  = 0.9,
  num_exemplars  = 300,
  do_plot        = FALSE
)

write.table(
  res_89$above_cutoff,
  file.path(out_dir, "similar_to_InsDel_N_beta_89.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

message("476-type: ", nrow(res_476$above_cutoff), " spectra above cutoff")
message("89-type: ",  nrow(res_89$above_cutoff),  " spectra above cutoff")
