# Rank signatures by average proportion of Del(T) from long poly-T
#
# Reads the partial credit table and prints signatures sorted by their
# average proportion across all tumors (descending).
#
# Usage:
#   Rscript rank_signatures.R

library(here)

pc <- read.delim(
  here::here(
    "code_for_internal_exploration/partial_credit_delT",
    "partial_credit_delT_long_homopolymers.tsv"
  ),
  check.names = FALSE
)

meta_cols <- c(
  "Patient", "Cancer_Type",
  "total_mutation_count", "del_T_long_poly_T_count", "del_T_long_poly_T_ratio"
)
sig_cols <- setdiff(colnames(pc), meta_cols)

sig_matrix <- as.matrix(pc[, sig_cols])
row_totals <- rowSums(sig_matrix)
row_totals[row_totals == 0] <- 1
sig_fractions <- sweep(sig_matrix, 1, row_totals, "/")

mean_prop <- sort(colMeans(sig_fractions), decreasing = TRUE)

result <- data.frame(
  Signature = names(mean_prop),
  Mean_Proportion = mean_prop,
  row.names = NULL
)

print(result, digits = 4, row.names = FALSE)
