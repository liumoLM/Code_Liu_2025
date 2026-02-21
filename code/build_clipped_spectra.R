#!/usr/bin/env Rscript
#
# Build clipped (R <= 9) spectra catalogs for all annotated VCF files
# from PCAWG and FMH datasets, for each of the three classification types
# (83, 89, 476). Writes one TSV per type to Manuscript_data/.

library(data.table)
library(ICAMS)
library(purrr)
library(stringr)

vcf_files <- sort(c(
  Sys.glob(path.expand(
    "~/MEGA/important_mut_sig_data/pcawg_indel_vcfs/*.annotated.indel.vcf.gz"
  )),
  Sys.glob(path.expand(
    "~/MEGA/important_mut_sig_data/fmh-unfiltered_vcfs/*.annotated.indel.vcf.gz"
  ))
))

stopifnot(length(vcf_files) > 0)
cat("Found", length(vcf_files), "VCF files\n")

sample_ids <- str_extract(basename(vcf_files), "^[^.]+")

catalog_fns <- list(
  "83" = ICAMS::annot_vcf_to_83_catalog,
  "89" = ICAMS::annot_vcf_to_89_catalog,
  "476" = ICAMS::annot_vcf_to_476_catalog
)

out_dir <- here::here("Manuscript_data")

for (type in names(catalog_fns)) {
  cat("Processing type", type, "...\n")
  fn <- catalog_fns[[type]]

  catalogs <- imap(set_names(vcf_files, sample_ids), \(path, sid) {
    vcf <- fread(path)
    message(path)
    if (nrow(vcf) == 0) {
      message("  Skipping ", sid, " (empty VCF)")
      return(NULL)
    }
    fn(vcf, sample_id = sid, clip_le_9 = TRUE)
  }) |>
    compact()

  combined <- do.call(cbind, catalogs)

  out_path <- file.path(out_dir, paste0(type, "_clipped_spectra.tsv"))
  write.table(combined, out_path, sep = "\t", quote = TRUE)
  cat("Wrote", out_path, "\n")
}
