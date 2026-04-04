library(data.table)
library(glue)

how_many_spectra = 2000

spectra <- read.delim(
  here::here("Manuscript_data/Liu_et_al_476_type_spectra.tsv"),
  check.names = FALSE,
  row.names = 1,
  sep = "\t"
)

# Identify non-::SP samples
non_sp <- !grepl("::SP", colnames(spectra))
spectra_fmh <- spectra[, non_sp]

# Top 200 by colSums
cs <- colSums(spectra_fmh)
top <- names(sort(cs, decreasing = TRUE))[1:how_many_spectra]
spectra_top <- spectra_fmh[, top]

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

spectra_top[del2_row, ] <- round(spectra_top[del2_row, ] * 0.1353)
spectra_top[delT_rows, ] <- round(spectra_top[delT_rows, ] * 0.2889 / 0.9355) # (0.2889 / 0.9355)) # 0.2889) <----
spectra_top[insT_rows, ] <- round(spectra_top[insT_rows, ] * 0.4420 / 0.9298)
spectra_top[delC_rows, ] <- round(spectra_top[delC_rows, ] * 0.4880 / 0.7093)
spectra_top[insC_rows, ] <- round(spectra_top[insT_rows, ] * 0.7180 / 0.6199)


out_path <- here::here("msi_study/fmh_converted_to_pcawg_spectra.tsv")
write.table(spectra_top, out_path, sep = "\t", quote = FALSE, col.names = NA)
cat(sprintf(
  "Saved %d rows x %d columns to %s\n",
  nrow(spectra_top),
  ncol(spectra_top),
  out_path
))
