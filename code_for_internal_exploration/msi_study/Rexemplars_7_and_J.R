library(gridExtra)
library(ggplot2)
library(philentropy)
library(mSigPlot)
library(glue)

source(here::here("code/find_many_similar.R"))

sigs = read.table(
  "Manuscript_data/Liu_et_al_final_89_type_signatures.tsv",
  sep = "\t",
  header = TRUE,
  row.names = 1
)


res1 <- find_many_similar(
  sig_path = here::here(
    "Manuscript_data/Liu_et_al_final_476_type_signatures.tsv"
  ),
  "InsDel_J",
  here::here("Manuscript_data/Liu_et_al_476_type_spectra.tsv"),
  cosine_cutoff = 0.9,
  num_exemplars = 300,
  out_pdf = "exemplars_for_IndDelJ_sim_0.9_num_300.pdf",
  min_mutations = 50
)
fwrite(res1$above_cutoff, here::here("msi_study/J_hits_ge_0.0.csv"))

res2 <- find_many_similar(
  sig_path = here::here(
    "Manuscript_data/Liu_et_al_final_476_type_signatures.tsv"
  ),
  "InsDel7",
  here::here("Manuscript_data/Liu_et_al_476_type_spectra.tsv"),
  cosine_cutoff = 0.9,
  num_exemplars = 300,
  out_pdf = "exemplars_for_IndDel7_sim_0.9_num_300.pdf",
  min_mutations = 50
)
fwrite(res2$above_cutoff, here::here("msi_study/7_hits_ge_0.0.csv"))
