source("code/plot_476_v3.R")
source("code/plot_476_pdf.R")
uu = read.delim("data/type476_spectra.tsv", sep = '\t')
zz = read.delim("data/type476_our_sigs.tsv", sep = '\t')


plot_476_pdf(uu[, c(3, 6, 19, 22, 100), drop = F], "foo1.pdf")

plot_476_pdf(zz[, c(3, 6, 19, 22, 100), drop = F], num_labels = 4, "foo2.pdf")
plot_476_pdf(
  zz[, c(3, 6, 19, 22, 100), drop = F],
  num_labels = 4,
  simplify_labels = TRUE,
  "foo3.pdf"
)
