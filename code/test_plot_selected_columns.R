source("code/plot_selected_columns.R")

# Plot specific signatures from an 83-channel file
plot_selected_columns(
  "data/type83_our_sigs.tsv",
  "tmp/selected_sigs_83.pdf",
  c("C_ID1$", "C_ID2$", "ID_N$")
)

dbs8 = read.csv("code/v3.2_DBS8_TISSUE_g7S68gD.txt", sep = '\t')

# Plot from 89-channel file
plot_selected_columns(
  "data/type89_spectra.tsv",
  "tmp/selected_spectra_89.pdf",
  dbs8$Sample.Names
)

# Plot from 476-channel file
plot_selected_columns(
  "data/type89_koh_sigs.tsv",
  "tmp/koh_89_all.pdf",
  ".*"
)

plot_selected_columns(
  "data/type89_our_sigs.tsv",
  "tmp/our_89_all.pdf",
  ".*"
)
