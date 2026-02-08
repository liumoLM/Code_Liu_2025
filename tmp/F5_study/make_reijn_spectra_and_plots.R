# make_reijn_spectra.R

library(magrittr)
library(ICAMS)
library(mSigPlot)
library(ggplot2)
library(cowplot)

# source("code/generate_89_and_476_type_catalog.R")
source("code/Generate_Koh89_Koh476_catalog_0121.R")

files = dir("Manuscript_data", pattern = "Reijn", full.names = TRUE)
mouse = files[1]
rpe1 = files[2]

mvcf = read.delim(mouse, sep = '\t')

rvcf = read.delim(rpe1, sep = '\t')

base_size83 <- 25
text_size83 <- 7
count_label_size83 <- 8

# MICE, type 476 and type 89
plot_mice = function() {
  m4 = GenerateKoh476CatalogfromAnnotateVcf(mvcf, "Sample")
  m_all = as.data.frame(rowSums(m4))
  p4_1 = plot_476(m_all, base_size = 20, plot_title = "Mutated mice")

  xm4 = GenerateKoh89CatalogfromAnnotateVcf(mvcf, "Sample")
  xm_all = as.data.frame(rowSums(xm4))
  p8_1 = plot_89(xm_all, base_size = 20, plot_title = "Mutated mice")

  colnames(mvcf)[1:2] <- c("CHROM", "POS")
  all_mice_83 = ICAMS::VCFsToIDCatalogs(
    list(all_mice = mvcf[, 1:8]),
    ref.genome = "mm10"
  )
  mouse83 = plot_83(
    all_mice_83[[1]],
    base_size = base_size83,
    text_size = text_size83,
    count_label_size = count_label_size83,
    plot_title = "Mutated mice"
  )
  return(c(mouse476 = p4_1, mouse89 = p8_1, mouse83 = mouse83))
}
mice = plot_mice()

plot_cells476 = function() {
  r4 = GenerateKoh476CatalogfromAnnotateVcf(rvcf, "Sample")
  r4_mut = r4[, 1:2]
  r4_wt = r4[, 3:5]
  r_all_mut = as.data.frame(rowSums(r4_mut))
  p4_2 = plot_476(r_all_mut, base_size = 20, plot_title = "Mutated cells")
  r_all_wt = as.data.frame(rowSums(r4_wt))
  p4_3 = plot_476(r_all_wt, base_size = 20, plot_title = "Wild-type cells")
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
  cell83_mut = plot_83(
    all_mut,
    plot_title = "Mutated cells",
    base_size = base_size83,
    text_size = text_size83,
    count_label_size = count_label_size83
  )
  all_wt = as.data.frame(rowSums(wt))
  cell83_wt = plot_83(
    all_wt,
    plot_title = "Wild-type cells",
    base_size = base_size83,
    text_size = text_size83,
    count_label_size = count_label_size83
  )
  return(c(cell83_mut, cell83_wt))
}
cells83 = plot_cells83()


cowplot::save_plot(
  "tmp/F5_study/insdelF_plots.pdf",
  gridExtra::marrangeGrob(
    grobs = c(
      mice["mouse476"],
      cells476,
      mice["mouse89"],
      cells89,
      mice["mouse83"],
      cells83
    ),
    ncol = 1,
    nrow = 6
  ),
  base_width = 20,
  base_height = 30
)
