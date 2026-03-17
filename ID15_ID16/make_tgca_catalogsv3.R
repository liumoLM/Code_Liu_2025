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

######

zzz = read.csv("TCGA-13-0889-01A-01W-0420-08_indels_categorized_vcf.csv")

make476catalog = function(tcgaid) {
  fullpath = dir(vcfdir, pattern = tcgaid, full.names = TRUE)
  b1 = ICAMS:::ReadVCF(fullpath, filter.status = NULL)
  b2 = ICAMS::AnnotateIDVCF(
    b1,
    ref.genome = "hg19",
    add_transcript_ranges = FALSE
  )
  b3 = cbind(b2$annotated.vcf, sample_id = tcgaid)
  return(GenerateKoh476CatalogfromAnnotateVcf(b3, "sample_id"))
}

foo = ICAMS::VCFsToCatalogs(
  "ID15_ID16/OVARIAN.SEROUS.CYSTADENOCARCINOMA.TCGA-13-0889-01A-01W-0420-08.GRCh37.indel.vcf",
  filter.status = "PASS",
  ref.genome = "hg19"
)
correct_id16 = foo$catID
correct_16_linking_tumor = foo$catID
colnames(correct_id16) = c('ID16')
correct_id16 <- correct_id16 / colSums(correct_id16)
write.table(
  correct_id16,
  "ID15_ID16/correct_ID16.tsv",
  sep = '\t',
  col.names = NA
)
uu = read.delim(
  "Manuscript_data/Liu_et_al_83_type_spectra.tsv",
  check.names = FALSE,
  row.names = 1
)
uu[, "WES::TCGA-13-0889-01A-01W-0420-08"] = correct_16_linking_tumor
write.table(
  uu,
  "Manuscript_data/Liu_et_al_83_type_spectra.tsv",
  sep = '\t',
  col.names = NA
)

bar = ICAMS::ReadVCFs(id16_vcf_path, filter.status = "PASS")
qq = ICAMS::AnnotateIDVCF(bar[[1]], ref.genome = "hg19")

# FROM HERE DOWN, IGNORE

# wget https://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz

idlist = c(
  "TCGA-EW-A2FV-01A-11D-A17D-09",
  "TCGA-D8-A27V-01A-12D-A17D-09",
  "TCGA-13-0889-01A-01W-0420-08"
)

cats476 = lapply(X = idlist, make476catalog)
names(cats476) = idlist
c476 = do.call(cbind, cats476)
write.table(c476, file = "ID15_I16_476.tsv", sep = '\t', row.names = TRUE)

cats89 = lapply(X = idlist, make89catalog)
names(cats89) = idlist
c89 = do.call(cbind, cats89)
write.table(c89, file = "ID15_I16_89.tsv", sep = '\t', row.names = TRUE)
colnames(c89)


source("~/github/Code_Liu_2025/code/Koh89_Koh476_Plotting_Functions.R")
library(PCAWG7)


plotall <- function(xx) {
  library(ggplotify)
  library(patchwork)

  c83 <- cats83[, xx, drop = FALSE]

  # Convert base R plots using as.grob()
  p1 <- as.grob(function() ICAMS::PlotCatalog(catalog = c83))

  p2 <- PlotKoh89Catalog(
    cats89[[xx]],
    text_size = 5,
    plot_title = "",
    setyaxis = NULL,
    ylabel = "Counts"
  )

  c476 <- cats476[[xx]]
  p3 <- PlotKoh476Catalog(
    c476,
    text_size = 5,
    plot_title = ""
  )

  s <- sbs[, xx, drop = FALSE]
  p4 <- as.grob(function() ICAMS::PlotCatalog(catalog = s))

  combined <- wrap_elements(p1) / p2 / p3 / wrap_elements(p4)

  filename <- file.path(pdfdir, paste0(xx, ".pdf"))
  ggsave(filename, plot = combined, width = 12, height = 16, device = "pdf")

  invisible(filename)
}

pdfdir = "~/github/Liu2024/misc-code-and-data"

lapply(idlist, plotall)

get_tgca_ = function(xxid) {
  xx = readr::read_delim(
    file.path(data_dir, "v0.2.7.PUBLIC.SORTED.DEDUP.UNIQ.simple"),
    delim = '\t',
    col_names = c(
      'cancer_type_long',
      'ID',
      'cancer_type',
      'genome_version',
      'mutation_type',
      'CHROM',
      'POS',
      'POS2',
      'REF',
      'ALT',
      'num_callers',
      'callers'
    ),
    col_type = "cccccciiccic"
  )
  # browser()
  # Filter to SNVs for SBS analysis
  xx = dplyr::filter(xx, ID == xxid, mutation_type == "SNV")

  # Create proper VCF-like data frame with required columns
  vcf_df = xx |>
    dplyr::select(CHROM, POS, REF, ALT) |>
    as.data.frame()

  yy = AnnotateSBSVCF(vcf_df, ref.genome = "hg19")
  write.table(yy, paste0(xxid, "_sbs.tsv"), sep = '\t')
  return(yy)
}

a1 = get_tgca_(idlist[1])
a2 = get_tgca_(idlist[2])
a1$id = idlist[1]
a2$id = idlist[2]
aa = rbind(a1, a2)
View(group_by(aa, trans.gene.symbol) %>% dplyr::filter(dplyr::n() > 1))


a3 = get_tgca_(idlist[3])


library(dplyr)
df <- read.csv("TCGA-13-0889-01A-01W-0420-08_sbs.csv", sep = "\t") # Adjust sep if needed

# Create the 5 required columns
annovar_input <- df %>%
  mutate(
    # Basic logic for SNVs; Indels may require more complex start/end adjustments
    Start = POS,
    End = POS
  ) %>%
  select(CHROM, Start, End, REF, ALT)
annovar_input$Otherinfo = "NONE"

write.table(
  annovar_input,
  "my_results.avinput",
  sep = "\t",
  quote = F,
  row.names = F,
  col.names = F
)

annovar_output = read.csv(
  "~/github/Code_Liu_2025/ID15_ID16/impact_results.hg19_multianno.txt",
  sep = '\t'
)

pqs = function(csv_vcf) {
  vcf = read.csv(csv_vcf)
  rr = list()
  # browser()
  for (i in 1:nrow(vcf)) {
    id = paste0("x_", vcf[i, ]$CHROM, vcf[i, ]$POS)
    sx = substring(vcf[i, ]$seq.context, 11, 110)
    px = pqsfinder::pqsfinder(DNAString(sx))
    rr[[id]] = pqsfinder::maxScores(px)
  }
  return(rr)
}

uu = pqs("TCGA-13-0889-01A-01W-0420-08_indels_categorized_vcf.csv")

create_tab_vcf = function(tcgaid) {
  zzz = read.csv(paste0(tcgaid, "_indels_categorized_vcf.csv"))
  zzz = zzz[, -1]
  colnames(zzz)[1] = '#CHROM'
  write.table(
    zzz,
    paste0(tcgaid, "_indels_categorized.vcf"),
    sep = "\t",
    row.names = FALSE
  )
}

lapply(idlist, create_tab_vcf)

plotall = function(tcgaid) {
  zzz = read.csv(paste0(tcgaid, "_indels_categorized_vcf.csv"))
  plot_pos_distances(zzz, bins = 10, max_dist = 1000)
}
lapply(idlist, plotall)
