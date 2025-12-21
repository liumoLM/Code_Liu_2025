# read_all_spectra

data_dir = "data"

type83_spectra = read.csv(
  file.path(data_dir, "type83_spectra.tsv"),
  sep = '\t'
)

type89_spectra = read.csv(
  file.path(data_dir, "type89_spectra.tsv"),
  sep = '\t'
)

type476_spectra = read.csv(
  file.path(data_dir, "type476_spectra.tsv"),
  sep = '\t'
)

stopifnot(colnames(type89_spectra) == colnames(type83_spectra))
stopifnot(colnames(type89_spectra) == colnames(type476_spectra))
