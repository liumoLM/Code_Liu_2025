# make_reijn_spectra.R

library(magrittr)
library(ICAMS)
library(mSigPlot)
library(ggplot2)
library(cowplot)

sigFfile = "DRUP01030028T"
sig4file = "CPCT02100089T"

# source("code/generate_89_and_476_type_catalog.R")
source("code/Generate_Koh89_Koh476_catalog_0121.R")

myplot83 = function(xx, plot_title) {
  plot_83(
    xx,
    plot_title = plot_title,
    base_size = base_size83,
    text_size = text_size83,
    count_label_size = count_label_size83
  )
}

myplot476 = function(xx, plot_title) {
  plot_476(
    xx,
    plot_title = plot_title,
    base_size = 8,
    label_size = 1
  )
}

files = dir("Manuscript_data", pattern = "Reijn", full.names = TRUE)
mouse = files[1]
rpe1 = files[2]

mvcf = read.delim(mouse, sep = '\t')

rvcf = read.delim(rpe1, sep = '\t')

base_size83 <- 25
text_size83 <- 7
count_label_size83 <- 6.5

# MICE, type 476 and type 89
plot_mice = function() {
  m4 = GenerateKoh476CatalogfromAnnotateVcf(mvcf, "Sample")
  m_all = as.data.frame(rowSums(m4))
  p4_1 = myplot476(m_all, plot_title = "Mutated mice")

  xm4 = GenerateKoh89CatalogfromAnnotateVcf(mvcf, "Sample")
  xm_all = as.data.frame(rowSums(xm4))
  p8_1 = plot_89(xm_all, base_size = 20, plot_title = "Mutated mice")

  colnames(mvcf)[1:2] <- c("CHROM", "POS")
  all_mice_83 = ICAMS::VCFsToIDCatalogs(
    list(all_mice = mvcf[, 1:8]),
    ref.genome = "mm10"
  )
  mouse83 = myplot83(all_mice_83[[1]], plot_title = "Mutated mice")
  return(c(mouse476 = p4_1, mouse89 = p8_1, mouse83 = mouse83))
}
mice = plot_mice()

plot_cells476 = function() {
  r4 = GenerateKoh476CatalogfromAnnotateVcf(rvcf, "Sample")
  r4_mut = r4[, 1:2]
  r4_wt = r4[, 3:5]
  r_all_mut = as.data.frame(rowSums(r4_mut))
  p4_2 = myplot476(r_all_mut, plot_title = "Mutated cells")
  r_all_wt = as.data.frame(rowSums(r4_wt))
  p4_3 = myplot476(r_all_wt, plot_title = "Wild-type cells")
  return(c(cell476_mut = p4_2, cell476_wt = p4_3))
}
cells476 = plot_cells476()

plot_cells89 = function() {
  xr4 = GenerateKoh89CatalogfromAnnotateVcf(rvcf, "Sample")
  xr4_mut = xr4[, 1:2]
  xr4_wt = xr4[, 3:5]
  xr_all_mut = as.data.frame(rowSums(xr4_mut))
  p8_2 = plot_89(xr_all_mut, base_size = 20, plot_title = "Mutated cells")
  xr_all_wt = as.data.frame(rowSums(xr4_wt))
  p8_3 = plot_89(xr_all_wt, base_size = 20, plot_title = "Wild-type cells")
  return(c(cell89_mut = p8_2, cell89_wt = p8_3))
}
cells89 = plot_cells89()

plot_cells83 = function() {
  rvcf %>%
    dplyr::select(CHROM, POS, ALT, REF, Sample) %>%
    dplyr::group_by(Sample) %>%
    dplyr::group_split() %>%
    ICAMS::VCFsToIDCatalogs(ref.genome = "hg38") %>%
    `$`(catalog) -> cats
  mut = cats[, 1:2]
  wt = cats[, 3:5]
  all_mut = as.data.frame(rowSums(mut))
  cell83_mut = myplot83(all_mut, plot_title = "Mutated cells")
  all_wt = as.data.frame(rowSums(wt))
  cell83_wt = myplot83(all_wt, plot_title = "Wild-type cells")
  return(c(cell83_mut, cell83_wt))
}
cells83 = plot_cells83()

## Plot linking tumors and signature spectra
stdread = function(ff) {
  read.delim(ff, sep = '\t', row.names = 1)
}

stdread("Manuscript_data/Liu_et_al_final_83_type_signatures.tsv") %>%
  `[`("ID_F") %>%
  myplot83("ID_F") -> sig83F

stdread("Manuscript_data/Liu_et_al_final_83_type_signatures.tsv") %>%
  `[`("C_ID4") %>%
  myplot83("C_ID4") -> sig834


# stdread("Manuscript_data/Liu_et_al_final_89_type_signatures.tsv") %>%
#  `[`("ID_F") %>%
#  plot_89() -> sigF_89

#stdread("Manuscript_data/Liu_et_al_final_89_type_signatures.tsv") %>%
#  `[`("C_ID4") %>%
#  plot_89() -> sig4_89

stdread("Manuscript_data/Liu_et_al_final_476_type_signatures.tsv") %>%
  dplyr::select(dplyr::contains("InsDel_F")) %>%
  myplot476("InsDelF") -> sigF476

stdread("Manuscript_data/Liu_et_al_final_476_type_signatures.tsv") %>%
  `[`("InsDel4") %>%
  myplot476("InsDel4") -> sig4476


stdread("Manuscript_data/Liu_et_al_476_type_spectra.tsv") %>%
  dplyr::select(dplyr::contains(sigFfile)) %>%
  myplot476(glue("InsDelF linking tumor {sigFfile}")) -> fF476

stdread("Manuscript_data/Liu_et_al_476_type_spectra.tsv") %>%
  dplyr::select(dplyr::contains(sig4file)) %>%
  myplot476(glue("InsDel4 linking tumor {sig4file}")) -> f4476

p_tight <- theme(plot.margin = margin(0, 0, 0, 0, "pt"))

grob1 <- gridExtra::arrangeGrob(
  grobs = c(
    sigF476 + p_tight,
    fF476 + p_tight,
    sig4476 + p_tight,
    f4476 + p_tight,
    mice["mouse476"] + p_tight,
    cells476["cell476_mut"] + p_tight,
    cells476["cell476_wt"] + p_tight
  ),
  ncol = 1,
  nrow = 7,
  padding = unit(0, "line")
)
ggsave(
  "F5_study/insdelF_plots.pdf",
  grob1,
  width = 10,
  height = 14,
  units = "in"
)

grob2 <- gridExtra::arrangeGrob(
  grobs = c(
    mice["mouse89"],
    cells89,
    mice["mouse83"],
    cells83,
    sig83F,
    sig834
  ),
  ncol = 1,
  padding = unit(0, "line")
)
ggsave(
  "F5_study/insdelF_more_plots.pdf",
  grob2,
  width = 10,
  height = 15,
  units = "in"
)
