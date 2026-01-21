# make_reijn_spectra.R

# source("code/generate_89_and_476_type_catalog.R")
source("code/Generate_Koh89_Koh476_catalog_0121.R")

files = dir("Manuscript_data", pattern = "Reijn", full.names = TRUE)
mouse = files[1]
rpe1 = files[2]

mvcf = read.delim(mouse, sep = '\t')

rvcf = read.delim(rpe1, sep = '\t')

library(mSigPlot)

mx = GenerateKoh476CatalogfromAnnotateVcf(mvcf, "Sample")
m4 = GenerateKoh476CatalogfromAnnotateVcf(mvcf, "Sample")
m_all = as.data.frame(rowSums(m4))
p4_1 = plot_476(m_all, base_size = 20, plot_title = "Mutated mice")

xm4 = GenerateKoh89CatalogfromAnnotateVcf(mvcf, "Sample")
xm_all = as.data.frame(rowSums(xm4))
p8_1 = plot_89(xm_all, base_size = 20, plot_title = "Mutated mice")

r4 = GenerateKoh476CatalogfromAnnotateVcf(rvcf, "Sample")
r4_mut = r4[, 1:2]
r4_wt = r4[, 3:5]
r_all_mut = as.data.frame(rowSums(r4_mut))
p4_2 = plot_476(r_all_mut, base_size = 20, plot_title = "Mutated cells")
r_all_wt = as.data.frame(rowSums(r4_wt))
p4_3 = plot_476(r_all_wt, base_size = 20, plot_title = "Wild-type cells")

xr4 = GenerateKoh89CatalogfromAnnotateVcf(rvcf, "Sample")
xr4_mut = xr4[, 1:2]
xr4_wt = xr4[, 3:5]
xr_all_mut = as.data.frame(rowSums(xr4_mut))
p8_2 = plot_89(xr_all_mut, base_size = 20, plot_title = "Mutated cells")
xr_all_wt = as.data.frame(rowSums(xr4_wt))
p8_3 = plot_89(xr_all_wt, base_size = 20, plot_title = "Wild-type cells")

ggsave(
  "tmp/F5_study/insdelF_plots.pdf",
  gridExtra::marrangeGrob(
    grobs = c(p4_1, p8_1, p4_2, p8_2, p4_3, p8_3),
    ncol = 1,
    nrow = 6
  ),
  width = 20,
  height = 30
)


plots = c()

for (cn in colnames(m4)) {
  p = plot_476(
    m4[, cn, drop = F],
    plot_title = paste("mouse", cn),
    base_size = 20
  )
  plots = c(plots, p)
}

ggsave(
  "mouse.pdf",
  gridExtra::marrangeGrob(grobs = plots, ncol = 1, nrow = 6),
  width = 20,
  height = 30
)


r4 = GenerateKoh476CatalogfromAnnotateVcf(rvcf, "Sample")
p2 = c()
for (cn in colnames(r4)) {
  p = plot_476(
    r4[, cn, drop = F],
    plot_title = paste("cell", cn),
    base_size = 20
  )
  p2 = c(p2, p)
}

ggsave(
  "cell.pdf",
  gridExtra::marrangeGrob(grobs = p2, ncol = 1, nrow = 6),
  width = 20,
  height = 30
)

mvcf %>%
  dplyr::filter(nchar(ins_or_del_seq) > 1) %>%
  dplyr::count(short_visual) %>%
  dplyr::arrange(desc(n)) %>%
  write.csv("reij_mouse.csv")

rvcf %>%
  dplyr::filter(nchar(ins_or_del_seq) > 1) %>%
  dplyr::count(short_visual) %>%
  dplyr::arrange(desc(n)) %>%
  write.csv("reij_rpe1.csv")
