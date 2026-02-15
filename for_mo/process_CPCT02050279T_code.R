source(here::here("code/Generate_Koh89_Koh476_catalog.R"))
########################################
###process purple somatic vcf###########
########################################

CPCT02050279T.somatic.vcf <- ICAMS::ReadVCFs(
  here::here("for_mo/CPCT02050279T.purple.somatic.vcf.gz"),
  variant.caller = "unknown",
  filter.status = "PASS"
)
this.vcf <- CPCT02050279T.somatic.vcf$CPCT02050279T.purple.somatic.vcf
this.vcf <- this.vcf[nchar(this.vcf$REF) != nchar(this.vcf$ALT), ] ##remove SNV and DBS, 204129 indels remaining


this.vcf.annotated <- ICAMS::AnnotateIDVCF(this.vcf, ref.genome = "hg19")

this.vcf.annotated.pass <- this.vcf.annotated$annotated.vcf #some complex indels were removed, 206101 indels remaining, some are duplicates because they were aligned to multiple transcripts

this.vcf.annotated.pass$mutID <- paste(
  "chr",
  this.vcf.annotated.pass$CHROM,
  "-",
  this.vcf.annotated.pass$POS,
  sep = ""
)
this.vcf.annotated.pass <- this.vcf.annotated.pass[
  !duplicated(this.vcf.annotated.pass$mutID),
] ## remove duplicates based on the chromosome and position of the indels, 203391 remaining

## edit the Koh476 annotation, change the long poly Ts to R(9,)
idx <- which(
  this.vcf.annotated.pass$R >= 9 &
    grepl(
      "Del\\(T\\)|Del\\(C\\)|Ins\\(C\\)|Ins\\(T\\)",
      this.vcf.annotated.pass$Koh_476
    )
)

strings <- this.vcf.annotated.pass$Koh_476[idx]
updated_strings <- gsub("R\\d+", "R(9,)", strings)
this.vcf.annotated.pass$Koh_476[idx] <- updated_strings

Koh476_catalog <- GenerateKoh476CatalogfromAnnotateVcf(
  this.vcf.annotated.pass,
  sample_col = "FILTER"
)

########################################
###process annotated somatic vcf###########
########################################

CPCT02050279T.annotated.vcf <- data.table::fread(
  here::here("for_mo/CPCT02050279T.annotated.indel.vcf.gz")
)
CPCT02050279T.annotated.vcf$mutID <- paste(
  "chr",
  CPCT02050279T.annotated.vcf$CHROM,
  "-",
  CPCT02050279T.annotated.vcf$POS,
  sep = ""
)
CPCT02050279T.annotated.vcf <- CPCT02050279T.annotated.vcf[
  CPCT02050279T.annotated.vcf$FILTER == "PASS",
] ## 206101 indels remaining, same as line 13

CPCT02050279T.annotated.vcf$mutID <- paste(
  "chr",
  CPCT02050279T.annotated.vcf$CHROM,
  "-",
  CPCT02050279T.annotated.vcf$POS,
  sep = ""
)
CPCT02050279T.annotated.vcf <- CPCT02050279T.annotated.vcf[
  !duplicated(CPCT02050279T.annotated.vcf$mutID),
] ## remove duplicates based on the chromosome and position of the indels, 203391 remaining

## edit the Koh476 annotation, change the long poly Ts to R(9,)
idx <- which(
  CPCT02050279T.annotated.vcf$R >= 9 &
    grepl(
      "Del\\(T\\)|Del\\(C\\)|Ins\\(C\\)|Ins\\(T\\)",
      CPCT02050279T.annotated.vcf$Koh_476
    )
)

strings <- CPCT02050279T.annotated.vcf$Koh_476[idx]
updated_strings <- gsub("R\\d+", "R(9,)", strings)
CPCT02050279T.annotated.vcf$Koh_476[idx] <- updated_strings

Koh476_catalog_from_annotated_indel_vcf <- GenerateKoh476CatalogfromAnnotateVcf(
  CPCT02050279T.annotated.vcf,
  sample_col = "FILTER"
)

sample_id <- "Colon::CPCT02050279T" # OK
spectra <- read.delim(
  here::here("Manuscript_data/Liu_et_al_476_type_spectra.tsv"),
  sep = '\t',
  check.names = FALSE,
  row.names = 1
)
sample_spectrum <- spectra[, sample_id, drop = FALSE]
