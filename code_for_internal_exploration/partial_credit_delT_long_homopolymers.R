# Partial credit of each signature to Del(T) from long T-homopolymers
#
# For each tumor, calculates each signature's contribution to deletions of
# 1 T from stretches of 6, 7, 8, or 9 T's using "partial credit":
#   partial_credit[sig, sample] = sum_over_target_channels(sig_476_profile[channel]) * assignment[sig, sample]
#
# Output:
#   /tmp/partial_credit_delT_long_homopolymers.tsv
#
# Usage:
#   Rscript partial_credit_delT_long_homopolymers.R

library(here)

data_dir <- here::here("Manuscript_data/finalized_cap9/")

# Load 89-type assignment matrix (signatures x samples)
assignments <- read.delim(
  file.path(data_dir, "liu_et_al_89_assignment.tsv"),
  row.names = 1,
  check.names = FALSE
)

# Load 476-type signature profiles (476 channels x signatures)
sigs_476 <- read.delim(
  file.path(data_dir, "liu_et_al_476_signatures.tsv"),
  row.names = 1,
  check.names = FALSE
)

# Load sample info for cancer type
sample_info <- read.delim(
  here::here("Manuscript_data", "sample_info.tsv"),
  check.names = FALSE
)

# Identify the 36 target channels: Del(T) from stretches of 6, 7, 8, or 9 T's
target_idx <- grep(
  "Del\\(T\\):R(6|7|8|\\(9,\\))",
  rownames(sigs_476)
)
target_channels <- rownames(sigs_476)[target_idx]
message("Found ", length(target_channels), " target channels")

# For each signature, sum its profile values across the 36 target channels
# This gives the fraction of each signature attributable to these Del(T) channels
target_fraction <- colSums(sigs_476[target_channels, , drop = FALSE])

# Align signatures between assignment matrix and signature profiles
common_sigs <- intersect(rownames(assignments), names(target_fraction))
message("Using ", length(common_sigs), " signatures common to both matrices")

assignments <- assignments[common_sigs, ]
target_fraction <- target_fraction[common_sigs]

# Compute partial credit matrix (sigs x samples):
# For each sig s and sample t: target_fraction[s] * assignment[s, t]
partial_credit <- sweep(
  as.matrix(assignments),
  1,
  target_fraction,
  "*"
)

# Transpose to samples x sigs, then build output data frame
result <- as.data.frame(t(partial_credit))
result$Patient <- rownames(result)

# Join with sample_info to get Cancer_Type
result <- merge(
  sample_info[, c("Patient", "Cancer_Type")],
  result,
  by = "Patient",
  all.y = TRUE
)
result$Cancer_Type[is.na(result$Cancer_Type)] <- "Unknown"

# Add summary columns
sig_cols <- common_sigs
result$total_mutation_count <- colSums(assignments)[result$Patient]
result$del_T_long_poly_T_count <- rowSums(result[, sig_cols])
result$del_T_long_poly_T_ratio <- result$del_T_long_poly_T_count /
  result$total_mutation_count

# Reorder: Patient, Cancer_Type, summary columns, then signatures
result <- result[, c(
  "Patient", "Cancer_Type",
  "total_mutation_count", "del_T_long_poly_T_count", "del_T_long_poly_T_ratio",
  sig_cols
)]

# Sort by Cancer_Type then Patient
result <- result[order(result$Cancer_Type, result$Patient), ]

out_path <- here::here(
  "code_for_internal_exploration/partial_credit_delT",
  "partial_credit_delT_long_homopolymers.tsv"
)
write.table(
  result,
  out_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
message("Wrote ", nrow(result), " samples to ", out_path)
