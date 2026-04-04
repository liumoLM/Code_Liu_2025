library(data.table)
library(glue)

j_hits <- read.csv(here::here("msi_study/J_hits_ge_0.9.csv"))
j_hits %>%
  dplyr::mutate(sampleid = gsub("..", "::", spectrum, fixed = TRUE)) -> jfixed

spectra <- read.delim(
  here::here("Manuscript_data/Liu_et_al_476_type_spectra.tsv"),
  check.names = FALSE,
  row.names = 1,
  sep = "\t"
)

spectra_top <- spectra[, jfixed$sampleid]

# Rows to scale
del2_row <- which(rownames(spectra_top) == "Del2:U1:R(5,9)")
delT_rows <- grep("Del\\(T\\):R\\(9,", rownames(spectra_top))
insT_rows <- grep("Ins\\(T\\):R\\(9", rownames(spectra_top))
insC_rows <- grep("Ins\\(C\\):R\\(9", rownames(spectra_top))
delC_rows <- grep("Del\\(C\\):R\\(9", rownames(spectra_top))

message(
  "Scaling\n",
  length(del2_row),
  " Del2:U1:R(5,9) rows\n",
  length(delT_rows),
  "Del(T):R(9,) rows, ",
  length(insT_rows),
  "Ins(T):R(9) rows\n",
  length(delC_rows),
  "Del(C):R(9,) rows, ",
  length(insC_rows),
  "Ins(C):R(9) rows\n"
)

#  From J_hits_ge_0.9.csv and 7_hits_ge_0.9.csv,
#  and then Rget_repeat_spectra.R which generated
#  check_[PCAWG_7|FMH_J]_Del2_U1_R5_9.csv
#
#  FMH_J (764,502 Del2:U1:R(5,9) mutations):
#  - R<=9: 57,078 (7.5%)
#  - R>9: 707,424 (92.5%)
#
#  PCAWG_7 (8,632 Del2:U1:R(5,9) mutations):
#  - R<=9: 6,097 (70.6%)
#  - R>9: 2,535 (29.4%)
old_spectra_top <- spectra_top
spectra_top[del2_row, ] <- round(
  spectra_top[del2_row, ] * 0.075 * (1 + 0.294 / 0.706)
)

spectra_top[delT_rows, ] <- round(spectra_top[delT_rows, ] * 0.3066 / 0.9435) # (0.2889 / 0.9355)) # 0.2889) <----
spectra_top[insT_rows, ] <- round(spectra_top[insT_rows, ] * 0.5032 / 0.9455)
spectra_top[delC_rows, ] <- round(spectra_top[delC_rows, ] * 0.5441 / 0.8049)
spectra_top[insC_rows, ] <- round(spectra_top[insT_rows, ] * 0.7738 / 0.9636)


out_path <- here::here("msi_study/fmh_J_converted_to_pcawg_spectra.tsv")
write.table(spectra_top, out_path, sep = "\t", quote = FALSE, col.names = NA)
cat(sprintf(
  "Saved %d rows x %d columns to %s\n",
  nrow(spectra_top),
  ncol(spectra_top),
  out_path
))
