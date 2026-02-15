vcf <- data.table::fread(
  "ID15_ID16/BREAST.INVASIVE.CARCINOMA.TCGA-D8-A27V-01A-12D-A17D-09.GRCh37.indel.vcf"
)
colnames(vcf)[1] <- "CHROM"

Rprof("aa.out")
annot_vcf <-
  ICAMS::AnnotateIDVCF(
    vcf[1:100, ],
    ref.genome = "hg19"
  )$annotated.vcf
Rprof(NULL)

source(here::here("code/annot_vcf_to_476_catalog.R"))

xx <- annot_vcf_to_476_catalog(annot_vcf)

newvcf <- data.table::fread(here::here(
  "for_mo/CPCT02050279T.annotated.indel.vcf.gz"
))

yy <- annot_vcf_to_476_catalog(newvcf)
