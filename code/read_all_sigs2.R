# read_all_sigs

data_dir = "data"

library(tidyr)


type89_our_sigs = read.delim(
  file.path(data_dir, "type89_liu_et_al_sigs.tsv"),
  sep = '\t'
)


type476_our_sigs = read.csv(
  file.path(data_dir, "type476_liu_et_al_sigs.tsv"),
  sep = '\t'
)

type89_koh_sigs = read.csv(
  file.path(data_dir, "type89_koh_sigs.tsv"),
  sep = '\t'
)

type83_our_sigs = read.delim(
  file.path(data_dir, "type83_liu_et_al_sigs.tsv"),
  sep = '\t'
)

type83_cosmic_sigs = read.delim(
  file.path(data_dir, "COSMIC_v3.5_ID_GRCh37.txt"),
  sep = '\t'
)

type83_jin_sigs = read.delim(
  file.path(data_dir, "jin_2024_indel_sigs_sup_tab_1.tsv"),
  sep = '\t'
)

all83 = cbind(type83_our_sigs, type83_cosmic_sigs, type83_jin_sigs)

all89 = cbind(type89_our_sigs, type89_koh_sigs)

linktable = read.delim(
  file.path(data_dir, "linktable.tsv"),
  sep = '\t'
)

library(dplyr)

us_v_all83 = read.csv(
  "signature_comparisons/us_vs_all_83.csv"
) |>
  select(-X, -max_cosine_id, -max_cosine, -euclidean.x, -euclidean.y)

us_v_koh89 = read.csv(
  "signature_comparisons/us_v_koh_89.csv"
) |>
  select(-euclidean)

link_and_83 = dplyr::full_join(
  linktable,
  us_v_all83,
  by = join_by(type83 == signature)
)

link_and_83_and_89 = dplyr::full_join(
  link_and_83,
  us_v_koh89,
  by = join_by(type89 == signature)
) |>
  rename(KohID = ID, cosine_v_koh = cosine)

us_v_spectra83 = read.delim("data/us_v_spectra_83.tsv", sep = '\t') |>
  select(-euclidean) |>
  rename(exemplar83 = ID, cosine_v_exemplar83 = cosine)

link_and_83_and_89_and_spectra83 = dplyr::full_join(
  link_and_83_and_89,
  us_v_spectra83,
  join_by(type83 == signature)
)

us_v_spectra89 = read.delim("data/us_v_spectra_89.tsv", sep = '\t') |>
  select(-euclidean) |>
  rename(exemplar89 = ID, cosine_v_exemplar89 = cosine)

link_and_83_and_89_and_spectra83_and_spectra89 = dplyr::full_join(
  link_and_83_and_89_and_spectra83,
  us_v_spectra89,
  join_by(type89 == signature)
)

View(link_and_83_and_89_and_spectra83_and_spectra89)

## To do, best_matches476; similar to best_matches89

## To do, plot everything
