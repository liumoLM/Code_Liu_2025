# make_reijn_spectra.R

source("code/generate_89_and_476_type_catalog.R")

files = dir("Manuscript_data", pattern = "Reijn", full.names = TRUE)
mouse = files[1]
rpe1 = files[2]

mvcf = read.delim(mouse, sep = '\t')

rvcf = read.delim(rpe1, sep = '\t')

library(mSigPlot)
library(Cairo)


m4 = GenerateKoh476CatalogfromAnnotateVcf(mvcf, "Sample")

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
