# Stacked bar charts of signature contributions to Del(T) from long poly-T
#
# For each cancer type, creates two sets of stacked bar charts:
#   Set 1: Top 30 tumors by del_T_long_poly_T_count
#   Set 2: Random sample of 30 from the remainder
# Each set has a counts plot (top) and proportions plot (bottom).
#
# Signatures contributing < min_contribution_fraction of any sample's total
# within the plotted subset are excluded from the legend.
#
# Arguments:
#   --min-contribution-fraction  Exclude signatures below this fraction (default: 0.02)
#   --max-tumors                 Max tumors per subset (default: 30)
#   --seed                       Random seed for reproducibility (default: 42)
#
# Output:
#   code_for_internal_exploration/partial_credit_delT/stacked_bars_delT.pdf
#
# Usage:
#   Rscript plot_stacked_bars.R
#   Rscript plot_stacked_bars.R --min-contribution-fraction 0.05

library(argparser)
library(ggplot2)
library(tidyr)
library(here)

p <- arg_parser("Stacked bar charts of Del(T) partial credit by cancer type")
p <- add_argument(p, "--min-contribution-fraction", type = "double", default = 0.02,
  help = "Exclude signatures contributing less than this fraction in a subset")
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

# Signature columns are everything after the 5 metadata columns
meta_cols <- c(
  "Patient", "Cancer_Type",
  "total_mutation_count", "del_T_long_poly_T_count", "del_T_long_poly_T_ratio"
)
sig_cols <- setdiff(colnames(pc), meta_cols)

out_path <- here::here(
  "code_for_internal_exploration/partial_credit_delT",
  "stacked_bars_delT.pdf"
)

cancer_types <- sort(unique(pc$Cancer_Type))
message("Generating plots for ", length(cancer_types), " cancer types")

#' Plot a counts + proportions stacked bar chart pair for a subset of tumors
#'
#' @param subset_data Data frame with Patient and signature columns, already
#'   ordered as desired for the x-axis.
#' @param sig_cols Character vector of signature column names.
#' @param title_suffix Character string appended to the cancer type in titles.
#' @param ct Character: cancer type name for the title.
#' @param min_frac Numeric: minimum fraction threshold for keeping a signature.
#' @return Invisible NULL. Draws to the current graphics device.
plot_pair <- function(subset_data, sig_cols, title_suffix, ct, min_frac) {
  if (nrow(subset_data) == 0 || all(rowSums(subset_data[, sig_cols]) == 0)) {
    return(invisible(NULL))
  }

  # Filter signatures by min contribution fraction within this subset
  sig_matrix <- as.matrix(subset_data[, sig_cols])
  row_totals <- rowSums(sig_matrix)
  row_totals[row_totals == 0] <- 1
  sig_fractions <- sweep(sig_matrix, 1, row_totals, "/")
  max_frac_per_sig <- apply(sig_fractions, 2, max)
  keep_sigs <- names(max_frac_per_sig)[max_frac_per_sig >= min_frac]

  if (length(keep_sigs) == 0) return(invisible(NULL))

  plot_df <- subset_data[, c("Patient", keep_sigs)]
  plot_df$Patient <- factor(plot_df$Patient, levels = subset_data$Patient)

  long_df <- tidyr::pivot_longer(
    plot_df,
    cols = -Patient,
    names_to = "Signature",
    values_to = "Count"
  )

  p_counts <- ggplot(long_df, aes(x = Patient, y = Count, fill = Signature)) +
    geom_bar(stat = "identity") +
    labs(
      title = paste0(ct, " — ", title_suffix, " (counts)"),
      x = NULL,
      y = "Mutation count"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
      legend.position = "right"
    )

  p_prop <- ggplot(long_df, aes(x = Patient, y = Count, fill = Signature)) +
    geom_bar(stat = "identity", position = "fill") +
    labs(
      title = paste0(ct, " — ", title_suffix, " (proportions)"),
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

  # Set 1: top N tumors
  top_data <- head(ct_data, args$max_tumors)
  plot_pair(top_data, sig_cols, "top tumors", ct, args$min_contribution_fraction)

  # Set 2: random sample of N from the remainder
  remainder <- ct_data[-(seq_len(min(args$max_tumors, nrow(ct_data)))), ]
  if (nrow(remainder) > 0) {
    n_sample <- min(args$max_tumors, nrow(remainder))
    random_idx <- sample(nrow(remainder), n_sample)
    random_data <- remainder[sort(random_idx), ]
    # Sort by descending count within this random subset
    random_data <- random_data[order(-random_data$del_T_long_poly_T_count), ]
    plot_pair(random_data, sig_cols, "random remaining tumors", ct,
              args$min_contribution_fraction)
  }
}

dev.off()
message("Wrote ", out_path)
