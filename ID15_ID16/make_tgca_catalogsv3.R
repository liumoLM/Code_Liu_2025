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

vvv = getannotvcf("TCGA-13-0889-01A-01W-0420-08")
vv = vvv$annotated.vcf

vv |>
  count(ins_or_del, short_visual) |>
  arrange(desc(n))

ICAMS::annot_vcf_to_83_catalog(vv)
