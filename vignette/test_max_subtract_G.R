# Self-contained test of max_subtract_signature for InsDel_G
# and its linking tumor DRUP01030015T

channel_names <- c(
  "[Del(C):R1]A", "[Del(C):R1]T", "[Del(C):R2]A", "[Del(C):R2]T",
  "[Del(C):R3]A", "[Del(C):R3]T", "[Del(C):R(4,5)]A", "[Del(C):R(4,5)]T",
  "[Del(C):R(1,5)]G", "Del(C):R(6,9)",
  "A[Del(T):R(1,4)]A", "A[Del(T):R(1,4)]C", "A[Del(T):R(1,4)]G",
  "C[Del(T):R(1,4)]A", "C[Del(T):R(1,4)]C", "C[Del(T):R(1,4)]G",
  "G[Del(T):R(1,4)]A", "G[Del(T):R(1,4)]C", "G[Del(T):R(1,4)]G",
  "A[Del(T):R(5,7)]A", "A[Del(T):R(5,7)]C", "A[Del(T):R(5,7)]G",
  "C[Del(T):R(5,7)]A", "C[Del(T):R(5,7)]C", "C[Del(T):R(5,7)]G",
  "G[Del(T):R(5,7)]A", "G[Del(T):R(5,7)]C", "G[Del(T):R(5,7)]G",
  "A[Del(T):R(8,)]A", "A[Del(T):R(8,)]C", "A[Del(T):R(8,)]G",
  "C[Del(T):R(8,)]A", "C[Del(T):R(8,)]C", "C[Del(T):R(8,)]G",
  "G[Del(T):R(8,)]A", "G[Del(T):R(8,)]C", "G[Del(T):R(8,)]G",
  "A[Ins(C):R0]A", "A[Ins(C):R0]T",
  "Ins(C):R(0,3)", "Ins(C):R(4,6)", "Ins(C):R(7,)",
  "A[Ins(T):R(0,4)]A", "A[Ins(T):R(0,4)]C", "A[Ins(T):R(0,4)]G",
  "C[Ins(T):R(0,4)]A", "C[Ins(T):R(0,4)]C", "C[Ins(T):R(0,4)]G",
  "G[Ins(T):R(0,4)]A", "G[Ins(T):R(0,4)]C", "G[Ins(T):R(0,4)]G",
  "A[Ins(T):R(5,7)]A", "A[Ins(T):R(5,7)]C", "A[Ins(T):R(5,7)]G",
  "C[Ins(T):R(5,7)]A", "C[Ins(T):R(5,7)]C", "C[Ins(T):R(5,7)]G",
  "G[Ins(T):R(5,7)]A", "G[Ins(T):R(5,7)]C", "G[Ins(T):R(5,7)]G",
  "A[Ins(T):R(8,)]A", "A[Ins(T):R(8,)]C", "A[Ins(T):R(8,)]G",
  "C[Ins(T):R(8,)]A", "C[Ins(T):R(8,)]C", "C[Ins(T):R(8,)]G",
  "G[Ins(T):R(8,)]A", "G[Ins(T):R(8,)]C", "G[Ins(T):R(8,)]G",
  "Del(2,4):R1", "Del(5,):R1",
  "Del(2,8):U(1,2):R(2,4)", "Del(2,):U(1,2):R(5,)",
  "Del(3,):U(3,):R2", "Del(3,):U(3,):R(3,)",
  "Ins(2,4):R0", "Ins(5,):R0", "Ins(2,4):R1", "Ins(5,):R1",
  "Ins(2,):R(2,4)", "Ins(2,):R(5,)",
  "Del(2,5):M1", "Del(3,5):M2", "Del(4,5):M(3,4)",
  "Del(6,):M1", "Del(6,):M2", "Del(6,):M3", "Del(6,):M(4,)",
  "Complex"
)

# InsDel_G signature (89-type, normalized proportions)
sig_to_subtract <- c(
  0.0182754980, 0.0396065450, 0.0157089500, 0.0147282720,
  0.0103379640, 0.0051837880, 0.0106656760, 0.0062765430,
  0.0072750390, 0.0073117230, 0.0049014280, 0.0065992280,
  0.0011249800, 0.0083849190, 0.0064231220, 0.0048626140,
  0.0028951190, 0.0057663980, 0.0042534460, 0.0031439750,
  0.0036193140, 0.0047254580, 0.0075946810, 0.0048103790,
  0.0043021650, 0.0029851300, 0.0031825370, 0.0019694040,
  0.0084129880, 0.0084772580, 0.0035020950, 0.0014469620,
  0.0056098620, 0.0018072540, 0.0013801100, 0.0038622720,
  0.0010121890, 0.0021590150, 0.0007123580, 0.0137420320,
  0.0082711940, 0.0036301320, 0.0053599050, 0.0059903010,
  0.0067723340, 0.0045483030, 0.0054518980, 0.0077778070,
  0.1049085470, 0.1122299340, 0.2006174460, 0.0098608100,
  0.0136210050, 0.0129391710, 0.0098743670, 0.0048554890,
  0.0054359910, 0.0111090530, 0.0084044680, 0.0138032280,
  0.0076585210, 0.0089616850, 0.0057921910, 0.0042629730,
  0.0027781100, 0.0028248470, 0.0047876510, 0.0043142010,
  0.0011535940, 0.0090245910, 0.0044794260, 0.0142794380,
  0.0015876540, 0.0081357910, 0.0029592850, 0.0091035310,
  0.0132761480, 0.0111078670, 0.0036615380, 0.0077101580,
  0.0068907800, 0.0056414250, 0.0022049020, 0.0011591240,
  0.0057161680, 0.0105454660, 0.0070062150, 0.0084766500,
  0.0000000000
)
names(sig_to_subtract) <- channel_names

# Linking tumor DRUP01030015T spectrum (89-type, raw counts)
spectrum <- c(
  14, 26, 12, 12, 8, 6, 8, 6, 11, 12,
  10, 5, 2, 6, 6, 5, 5, 2, 5, 1,
  0, 6, 6, 1, 5, 2, 5, 0, 6, 9,
  2, 3, 6, 2, 3, 1, 1, 0, 2, 9,
  6, 5, 6, 3, 4, 4, 5, 5, 79, 95,
  128, 14, 13, 20, 11, 8, 7, 11, 9, 9,
  20, 19, 19, 12, 13, 13, 6, 13, 10, 19,
  15, 18, 6, 9, 11, 20, 25, 9, 10, 6,
  24, 17, 10, 4, 8, 10, 8, 6, 0
)
names(spectrum) <- channel_names

# Run the search loop as in vhelpers.R
min_prob_ge_total_negative <- 0.7

cat("=== Testing max_subtract_signature for InsDel_G ===\n")
cat("Spectrum: DRUP01030015T (total mutations:",
    sum(spectrum), ")\n\n")

best_max_neg_fraction <- NULL
best_mss_result <- NULL

for (max_neg_fraction in seq(0.005, 0.15, by = 0.005)) {
  mss_result <- mSigBG::max_subtract_signature(
    spectrum = spectrum,
    sig_to_subtract = sig_to_subtract,
    max_neg_fraction = max_neg_fraction
  )
  cat(sprintf(
    paste0(
      "max_neg_fraction = %.3f  =>  n_subtract = %7.1f, ",
      "n_residual = %7.1f, total_negative = %7.1f, ",
      "n_neg_channels = %2d, prob_ge_total_neg = %.3f%s\n"
    ),
    max_neg_fraction,
    mss_result$n_subtract,
    mss_result$n_residual,
    mss_result$total_negative,
    mss_result$n_negative_channels,
    mss_result$prob_ge_total_negative,
    if (mss_result$prob_ge_total_negative >= min_prob_ge_total_negative)
      "  *" else ""
  ))
  if (mss_result$prob_ge_total_negative >= min_prob_ge_total_negative) {
    best_max_neg_fraction <- max_neg_fraction
    best_mss_result <- mss_result
  }
}

cat("\n=== Result ===\n")
if (is.null(best_mss_result)) {
  cat("No max_neg_fraction met the threshold.\n")
  cat("Falling back to max_neg_fraction = 0.005\n")
  best_mss_result <- mSigBG::max_subtract_signature(
    spectrum = spectrum,
    sig_to_subtract = sig_to_subtract,
    max_neg_fraction = 0.005
  )
  best_max_neg_fraction <- 0.005
}

cat(sprintf("Best max_neg_fraction: %.3f\n", best_max_neg_fraction))
cat(sprintf("N subtract: %.1f\n", best_mss_result$n_subtract))
cat(sprintf("N residual: %.1f\n", best_mss_result$n_residual))
cat(sprintf("Total negative: %.1f\n", best_mss_result$total_negative))
cat(sprintf("N negative channels: %d\n", best_mss_result$n_negative_channels))
cat(sprintf("P(>= total negative): %.3f\n",
            best_mss_result$prob_ge_total_negative))

# Compute partial spectrum and residual
partial_spectrum <- best_mss_result$n_subtract * sig_to_subtract
residual_spectrum <- spectrum - partial_spectrum

# Format as named single-column matrices for mSigPlot::plot_89
as_catalog <- function(vec, col_name) {
  m <- matrix(vec, ncol = 1, dimnames = list(channel_names, col_name))
  m
}

partial_cat <- as_catalog(partial_spectrum, "InsDel_G partial spectrum")
residual_cat <- as_catalog(residual_spectrum, "Residual")

# Use a common y-axis so the two plots are visually comparable
ymax <- max(c(partial_spectrum, residual_spectrum))

library(mSigPlot)
library(ggplot2)

p_partial <- plot_89(
  partial_cat,
  plot_title = sprintf(
    "InsDel_G partial spectrum (n_subtract = %.1f, max_neg_fraction = %.3f)",
    best_mss_result$n_subtract, best_max_neg_fraction
  ),
  setyaxis = ymax
)

p_residual <- plot_89(
  residual_cat,
  plot_title = sprintf(
    "Residual (n_residual = %.1f, total_negative = %.1f, P >= total neg = %.3f)",
    best_mss_result$n_residual, best_mss_result$total_negative,
    best_mss_result$prob_ge_total_negative
  ),
  setyaxis = ymax
)

# Save plots to /tmp and open them
ggsave("/tmp/InsDel_G_partial_spectrum.png", p_partial,
       width = 10, height = 4, dpi = 150, bg = "white")
ggsave("/tmp/InsDel_G_residual.png", p_residual,
       width = 10, height = 4, dpi = 150, bg = "white")

cat("\nPlots saved to:\n")
cat("  /tmp/InsDel_G_partial_spectrum.png\n")
cat("  /tmp/InsDel_G_residual.png\n")
