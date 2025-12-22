# read_all_sigs

data_dir = "data"

library(tidyr)

s11 = read.csv(
  file.path(data_dir, "s11_type83.csv")
)

s11 %>%
  unite(
    "ID",
    Type,
    Subtype,
    Indel_size,
    Repeat_MH_size,
    sep = ":"
  ) -> xx

rownames(xx) = xx$ID
write.table(xx[, -1], file.path("data", "type83_our_sigs.tsv"), sep = '\t')

type83_our_sigs = read.csv(
  file.path(data_dir, "type83_our_sigs.tsv"),
  sep = '\t'
)
all(rownames(type83_spectra) == rownames(type83_our_sigs))
rm(xx, s11)

s12 = read.csv(
  file.path(data_dir, "s12_type89.csv")
)
rownames(s12) = s12[, 1]
write.table(s12[, -1], file.path(data_dir, "type89_our_sigs.tsv"), sep = '\t')

type89_our_sigs = read.csv(
  file.path(data_dir, "type89_our_sigs.tsv"),
  sep = '\t'
)
all(rownames(type89_our_sigs) == rownames(type89_spectra))
rm(s12)


s13 = read.csv(
  file.path(data_dir, "s13_type476.csv")
)
rownames(s13) = s13[, 1]
write.table(s13[, -1], file.path(data_dir, "type476_our_sigs.tsv"), sep = '\t')

type476_our_sigs = read.csv(
  file.path(data_dir, "type476_our_sigs.tsv"),
  sep = '\t'
)
all(rownames(type476_our_sigs) == rownames(type476_spectra))
rm(s13)

koh_s10 = read.csv(
  file.path(data_dir, "koh_s10_consensus_sigs.csv"),
)
rownames(koh_s10) = koh_s10[, 1]
write.table(
  koh_s10[, -1],
  file.path(data_dir, "type89_koh_sigs.tsv"),
  sep = '\t'
)
type89_koh_sigs = read.csv(
  file.path(data_dir, "type89_koh_sigs.tsv"),
  sep = '\t'
)
