# One-off script to create indel catalogs from the "simple" file

library(ICAMS)
library(readr)
library(tidyr)
library(stringi)
library(dplyr)
library(ICAMS)
data_dir = "~/MEGA/important_mut_sig_data/tcga_from_2020_paper"

vcfdir = file.path(data_dir, "vcfs")

getannotvcf = function(tcgaid) {
  fullpath = dir(vcfdir, pattern = tcgaid, full.names = TRUE)
  b1 = ICAMS:::ReadVCF(fullpath, filter.status = NULL)
  b2 = ICAMS::AnnotateIDVCF(
    b1,
    ref.genome = "hg19",
    add_transcript_ranges = FALSE,
    context_width_multiplier = 60L
  )

  return(b2)
}

tcgaid <- "TCGA-13-0889-01A-01W-0420-08"
vvv = getannotvcf(tcgaid)
vv = vvv$annotated.vcf

vv |>
  count(ins_or_del, short_visual) |>
  arrange(desc(n))

t83 <- ICAMS::annot_vcf_to_83_catalog(vv, sample_id = tcgaid)
t89 <- ICAMS::annot_vcf_to_89_catalog(vv, sample_id = tcgaid)
t476 <- ICAMS::annot_vcf_to_476_catalog(vv, sample_id = tcgaid)

mydir <- here::here("ID15_ID16/correct_id16_spectra")
write.table(t83, file.path(mydir, glue::glue("{tcgaid}_83.tsv")), sep = '\t')
write.table(t89, file.path(mydir, glue::glue("{tcgaid}_89.tsv")), sep = '\t')
write.table(t476, file.path(mydir, glue::glue("{tcgaid}_476.tsv")), sep = '\t')
