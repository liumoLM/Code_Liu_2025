# Demonstration of assignment vs spectra mismatch
# The 89-type assignment totals are exactly 2x the 89-type spectra totals.

library(here)
library(mSigPlot)
library(lsa)

source(here::here("vignette", "vhelpers.R"))

spectra_89 <- read_finalized("89_spectra")
assign_89 <- read_finalized("89_assignment")
sigs_89 <- read_finalized("89_signatures")

# ============================================================
# Problem 1: Assignment totals are 2x the spectra totals
# ============================================================

# First 5 shared columns
shared <- intersect(colnames(spectra_89), colnames(assign_89))
cols5 <- shared[1:5]

cat("=== Column names match ===\n")
cat("Spectra columns: ", paste(cols5, collapse = ", "), "\n")
cat("Assignment columns:", paste(cols5, collapse = ", "), "\n\n")

cat("=== colSums comparison (first 5 samples) ===\n")
spec_sums <- colSums(spectra_89[, cols5])
asgn_sums <- colSums(assign_89[, cols5])
comparison <- data.frame(
  sample = cols5,
  spectra_total = spec_sums,
  assign_total = asgn_sums,
  ratio = asgn_sums / spec_sums
)
print(comparison, row.names = FALSE)

# ============================================================
# Problem 2: InsDel1c — exemplar SP112907 has 0% InsDel1c
# ============================================================

cat("\n\n=== InsDel1c: exemplar SP112907 ===\n")

exemplar_1c <- "SP112907"

# Cosine similarity
cos_1c <- cosine(
  as.numeric(sigs_89[, "InsDel1c"]),
  as.numeric(spectra_89[, exemplar_1c])
)
cat(
  "Cosine similarity (InsDel1c sig vs",
  exemplar_1c,
  "spectrum):",
  cos_1c,
  "\n\n"
)

# Assignment for this exemplar
cat("Assignments for", exemplar_1c, ":\n")
a_1c <- assign_89[, exemplar_1c]
names(a_1c) <- rownames(assign_89)
print(a_1c[a_1c > 0])
cat("Total assigned:", sum(a_1c), "\n")
cat(
  "InsDel1c assigned:",
  a_1c["InsDel1c"],
  sprintf("(%.1f%%)\n", 100 * a_1c["InsDel1c"] / sum(a_1c))
)

# Plot signature and exemplar spectrum
sig_plot <- mSigPlot::plot_89(
  sigs_89[, "InsDel1c", drop = FALSE],
  plot_title = "InsDel1c signature"
)
spec_plot <- mSigPlot::plot_89(
  spectra_89[, exemplar_1c, drop = FALSE],
  plot_title = paste0(
    exemplar_1c,
    " spectrum (cosine = ",
    format(cos_1c, digits = 4),
    ")"
  )
)
gridExtra::grid.arrange(sig_plot, spec_plot, ncol = 1)

# ============================================================
# Problem 3: InsDel2a — exemplar CPCT02010955T is only 52%
# ============================================================

cat("\n\n=== InsDel2a: exemplar CPCT02010955T ===\n")

exemplar_2a <- "CPCT02010955T"

# Cosine similarity
cos_2a <- cosine(
  as.numeric(sigs_89[, "InsDel2a"]),
  as.numeric(spectra_89[, exemplar_2a])
)
cat(
  "Cosine similarity (InsDel2a sig vs",
  exemplar_2a,
  "spectrum):",
  cos_2a,
  "\n\n"
)

# Assignment for this exemplar
cat("Assignments for", exemplar_2a, ":\n")
a_2a <- assign_89[, exemplar_2a]
names(a_2a) <- rownames(assign_89)
print(a_2a[a_2a > 0])
cat("Total assigned:", sum(a_2a), "\n")
cat("Spectra total:", sum(spectra_89[, exemplar_2a]), "\n")
cat(
  "InsDel2a assigned:",
  a_2a["InsDel2a"],
  sprintf("(%.1f%%)\n", 100 * a_2a["InsDel2a"] / sum(a_2a))
)

# Plot signature and exemplar spectrum

sig_plot <- mSigPlot::plot_89(
  sigs_89[, "InsDel2a", drop = FALSE],
  plot_title = "InsDel2a signature"
)
spec_plot <- mSigPlot::plot_89(
  spectra_89[, exemplar_2a, drop = FALSE],
  plot_title = paste0(
    exemplar_2a,
    " spectrum (cosine = ",
    format(cos_2a, digits = 4),
    ")"
  )
)
gridExtra::grid.arrange(sig_plot, spec_plot, ncol = 1)

# ============================================================
# Problem 4: InsDel_B — exemplar CPCT02020413T
# ============================================================

cat("\n\n=== InsDel_B: exemplar CPCT02020413T ===\n")

exemplar_B <- "CPCT02020413T"
tsig = "InsDel_B"

# Cosine similarity
cos_B <- cosine(
  as.numeric(sigs_89[, tsig]),
  as.numeric(spectra_89[, exemplar_B])
)
cat(
  "Cosine similarity (InsDel_B sig vs",
  exemplar_B,
  "spectrum):",
  cos_B,
  "\n\n"
)

# Assignment for this exemplar
cat("Assignments for", exemplar_B, ":\n")
a_b <- assign_89[, exemplar_2a]
names(a_b) <- rownames(assign_89)
print(a_b[a_b > 0])
cat("Total assigned:", sum(a_b), "\n")
cat("Spectra total:", sum(spectra_89[, exemplar_B]), "\n")
cat(
  "InsDel_B assigned:",
  a_b[tsig],
  sprintf("(%.1f%%)\n", 100 * a_b[tsig] / sum(a_b))
)

# Plot signature and exemplar spectrum

sig_plot <- mSigPlot::plot_89(
  sigs_89[, tsig, drop = FALSE],
  plot_title = paste(tsig, "signature")
)
spec_plot <- mSigPlot::plot_89(
  spectra_89[, exemplar_B, drop = FALSE],
  plot_title = paste0(
    exemplar_B,
    " spectrum (cosine = ",
    format(cos_B, digits = 4),
    ")"
  )
)
gridExtra::grid.arrange(sig_plot, spec_plot, ncol = 1)
