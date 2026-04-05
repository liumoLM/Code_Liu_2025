# Stacked bar charts of signature contributions to Del(T) from long poly-T
#
# For each cancer type, creates two sets of stacked bar charts:
#   Set 1: Top 30 tumors by del_T_long_poly_T_count
#   Set 2: Random sample of 30 from the remainder
# Each set has a counts plot (top) and proportions plot (bottom).
#
# Signatures are grouped into 3 tiers of 6 by global average proportion:
#   Tier 1 (top 6):  bold saturated colors
#   Tier 2 (next 6): medium-tone colors
#   Tier 3 (next 6): pastel/light colors
#   Remaining:       grey "Other"
# Legend labels include tier prefix (T1/T2/T3) for quick identification.
#
# Arguments:
#   --max-tumors  Max tumors per subset (default: 30)
#   --seed        Random seed for reproducibility (default: 42)
#
# Output:
#   code_for_internal_exploration/partial_credit_delT/stacked_bars_delT.pdf
#
# Usage:
#   Rscript plot_stacked_bars.R

library(argparser)
library(ggplot2)
library(tidyr)
library(here)

p <- arg_parser("Stacked bar charts of Del(T) partial credit by cancer type")
p <- add_argument(p, "--max-tumors", type = "integer", default = 30L,
  help = "Maximum number of tumors per subset")
p <- add_argument(p, "--seed", type = "integer", default = 42L,
  help = "Random seed for the remainder sample")
args <- parse_args(p)

set.seed(args$seed)

# Read partial credit table
pc <- read.delim(
  here::here(
    "code_for_internal_exploration/partial_credit_delT",
    "partial_credit_delT_long_homopolymers.tsv"
  ),
  check.names = FALSE
)

# Join MSI status from sample_info
sample_info <- read.delim(
  here::here("Manuscript_data", "sample_info.tsv"),
  check.names = FALSE
)
pc$MSIseq_MSI.H <- sample_info$MSIseq_MSI.H[
  match(pc$Patient, sample_info$Patient)
]

meta_cols <- c(
  "Patient", "Cancer_Type",
  "total_mutation_count", "del_T_long_poly_T_count", "del_T_long_poly_T_ratio",
  "MSIseq_MSI.H"
)
sig_cols <- setdiff(colnames(pc), meta_cols)

# Rank signatures by global average proportion
sig_matrix <- as.matrix(pc[, sig_cols])
row_totals <- rowSums(sig_matrix)
row_totals[row_totals == 0] <- 1
sig_fractions <- sweep(sig_matrix, 1, row_totals, "/")
mean_prop <- sort(colMeans(sig_fractions), decreasing = TRUE)
ranked_sigs <- names(mean_prop)

# Define tiers
tier1_sigs <- ranked_sigs[1:6]   # bold saturated
tier2_sigs <- ranked_sigs[7:12]  # medium tone
tier3_sigs <- ranked_sigs[13:18] # pastel/light
other_sigs <- ranked_sigs[19:length(ranked_sigs)]

# Three visually distinct color families (6 each)
# Tier 1: bold saturated
tier1_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628")
# Tier 2: medium tones (distinct from tier 1)
tier2_colors <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F")
# Tier 3: pastels/light (clearly lighter)
tier3_colors <- c("#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#E5D8BD")

# Stacking order: Other at bottom, tier3 above, tier2 above, tier1 on top
all_levels <- c("Other", rev(tier3_sigs), rev(tier2_sigs), rev(tier1_sigs))

# Build fill color map with tier-prefixed legend labels
fill_map <- c(
  setNames("grey70", "Other"),
  setNames(rev(tier3_colors), rev(tier3_sigs)),
  setNames(rev(tier2_colors), rev(tier2_sigs)),
  setNames(rev(tier1_colors), rev(tier1_sigs))
)

out_path <- here::here(
  "code_for_internal_exploration/partial_credit_delT",
  "stacked_bars_delT.pdf"
)

cancer_types <- sort(unique(pc$Cancer_Type))
message("Generating plots for ", length(cancer_types), " cancer types")
message("Tier 1 (bold):    ", paste(tier1_sigs, collapse = ", "))
message("Tier 2 (medium):  ", paste(tier2_sigs, collapse = ", "))
message("Tier 3 (pastel):  ", paste(tier3_sigs, collapse = ", "))
message("Other: ", length(other_sigs), " signatures lumped")

plot_pair <- function(subset_data, title_suffix, ct,
                      tier1_sigs, tier2_sigs, tier3_sigs, other_sigs,
                      all_levels, fill_map) {
  if (nrow(subset_data) == 0) return(invisible(NULL))

  named_sigs <- c(tier1_sigs, tier2_sigs, tier3_sigs)
  plot_df <- subset_data[, c("Patient", named_sigs), drop = FALSE]
  plot_df$Other <- rowSums(subset_data[, other_sigs, drop = FALSE])
  plot_df$Patient <- factor(plot_df$Patient, levels = subset_data$Patient)

  value_cols <- c(named_sigs, "Other")

  long_df <- tidyr::pivot_longer(
    plot_df,
    cols = all_of(value_cols),
    names_to = "Signature",
    values_to = "Count"
  )
  long_df$Signature <- factor(long_df$Signature, levels = all_levels)

  # Only show signatures with nonzero contribution in the legend
  present_sigs <- names(which(colSums(plot_df[, value_cols, drop = FALSE]) > 0))
  legend_breaks <- intersect(all_levels, present_sigs)

  p_counts <- ggplot(long_df, aes(x = Patient, y = Count, fill = Signature)) +
    geom_col(colour = "grey30", linewidth = 0.15) +
    scale_fill_manual(values = fill_map, breaks = legend_breaks) +
    labs(
      title = paste0(ct, " \u2014 ", title_suffix, " (counts)"),
      x = NULL,
      y = "Mutation count"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
      legend.position = "right"
    )

  p_prop <- ggplot(long_df, aes(x = Patient, y = Count, fill = Signature)) +
    geom_col(position = "fill", colour = "grey30", linewidth = 0.15) +
    scale_fill_manual(values = fill_map, breaks = legend_breaks) +
    labs(
      title = paste0(ct, " \u2014 ", title_suffix, " (proportions)"),
      x = NULL,
      y = "Proportion"
    ) +
    scale_y_continuous(labels = scales::percent) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
      legend.position = "right"
    )

  gridExtra::grid.arrange(p_counts, p_prop, ncol = 1)
}

pdf(out_path, width = 14, height = 10)

for (ct in cancer_types) {
  ct_data <- pc[pc$Cancer_Type == ct, ]
  ct_data <- ct_data[order(-ct_data$del_T_long_poly_T_count), ]

  # Set 1: MSI-H tumors
  msi_data <- ct_data[!is.na(ct_data$MSIseq_MSI.H) & ct_data$MSIseq_MSI.H == TRUE, ]
  if (nrow(msi_data) > 0) {
    plot_pair(msi_data, paste0("MSI-H (n=", nrow(msi_data), ")"), ct,
              tier1_sigs, tier2_sigs, tier3_sigs, other_sigs,
              all_levels, fill_map)
  }

  # Non-MSI-H tumors
  non_msi <- ct_data[is.na(ct_data$MSIseq_MSI.H) | ct_data$MSIseq_MSI.H == FALSE, ]
  non_msi <- non_msi[order(-non_msi$del_T_long_poly_T_count), ]

  # Set 2: top N non-MSI-H tumors
  if (nrow(non_msi) > 0) {
    top_non_msi <- head(non_msi, args$max_tumors)
    plot_pair(top_non_msi, "top non-MSI-H", ct,
              tier1_sigs, tier2_sigs, tier3_sigs, other_sigs,
              all_levels, fill_map)
  }

  # Set 3: random N from remaining non-MSI-H
  remainder <- non_msi[-(seq_len(min(args$max_tumors, nrow(non_msi)))), ]
  if (nrow(remainder) > 0) {
    n_sample <- min(args$max_tumors, nrow(remainder))
    random_idx <- sample(nrow(remainder), n_sample)
    random_data <- remainder[sort(random_idx), ]
    random_data <- random_data[order(-random_data$del_T_long_poly_T_count), ]
    plot_pair(random_data, "random non-MSI-H", ct,
              tier1_sigs, tier2_sigs, tier3_sigs, other_sigs,
              all_levels, fill_map)
  }
}

dev.off()
message("Wrote ", out_path)
