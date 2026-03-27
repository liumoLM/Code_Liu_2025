# Demonstration of assignment vs spectra mismatch
# The 89-type assignment totals are exactly 2x the 89-type spectra totals.

library(here)
library(mSigPlot)
library(lsa)

source(here::here("vignette", "vhelpers.R"))

spectra_89 <- read_finalized("89_spectra")
assign_89 <- read_finalized("89_assignment")
sigs_89 <- read_finalized("89_signatures")
sig <- "InsDel1b"

exemplar <- "SP21400"

# Cosine similarity
cos_1c <- cosine(
  as.numeric(sigs_89[, sig]),
  as.numeric(spectra_89[, exemplar])
)
cat(
  "Cosine similarity (InsDel1c sig vs",
  exemplar,
  "spectrum):",
  cos_1c,
  "\n\n"
)

# Assignment for this exemplar
cat("Assignments for", exemplar, ":\n")
a_1c <- assign_89[, exemplar]
names(a_1c) <- rownames(assign_89)
print(a_1c[a_1c > 0])
cat("Total assigned:", sum(a_1c), "\n")
cat(
  sig,
  "assigned:",
  a_1c[sig],
  sprintf("(%.1f%%)\n", 100 * a_1c[sig] / sum(a_1c))
)

# Plot signature and exemplar spectrum
sig_plot <- mSigPlot::plot_89(
  sigs_89[, sig, drop = FALSE],
  plot_title = sig
)
spec_plot <- mSigPlot::plot_89(
  spectra_89[, exemplar, drop = FALSE],
  plot_title = paste0(
    exemplar,
    " spectrum (cosine = ",
    format(cos_1c, digits = 4),
    ")"
  )
)
gridExtra::grid.arrange(sig_plot, spec_plot, ncol = 1)
