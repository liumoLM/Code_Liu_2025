## Triangular correlation heatmaps for signature co-occurrence
##
## Input: liu_et_al_89_assignment.tsv (signatures × samples)
##        sample_info.tsv (sample metadata with Cancer_Type)
## Output: One multi-page PDF with:
##   1. Spearman correlation on log(absolute counts + 1) - all samples
##   2. Pearson correlation on signature proportions - all samples
##   3-N. Same two heatmaps per cancer type

library(ggplot2)
library(reshape2)

assignment_file <- here::here(
  "Manuscript_data/finalized_cap9/liu_et_al_89_assignment.tsv"
)
sample_info_file <- here::here("Manuscript_data/sample_info.tsv")

df <- read.delim(assignment_file, row.names = 1, check.names = FALSE)
sample_info <- read.delim(sample_info_file, check.names = FALSE)

# Map sample IDs to cancer types
cancer_map <- setNames(sample_info$Cancer_Type, sample_info$Patient)

# Helper: compute both correlation matrices from a signature × sample matrix
compute_cors <- function(mat) {
  # Spearman on log(counts + 1)
  log_counts <- log1p(as.matrix(mat))
  cor_sp <- cor(t(log_counts), method = "spearman")

  # Pearson on proportions
  cs <- colSums(mat)
  props <- sweep(as.matrix(mat), 2, cs, "/")
  props <- props[, is.finite(colSums(props)), drop = FALSE]
  cor_pe <- cor(t(props), method = "pearson")

  list(spearman = cor_sp, pearson = cor_pe)
}

# Helper to make a lower-triangle heatmap
make_tri_heatmap <- function(cor_mat, title) {
  sig_names <- rownames(cor_mat)
  # Keep strictly lower triangle (exclude diagonal)
  cor_mat[upper.tri(cor_mat, diag = TRUE)] <- NA

  melted <- melt(cor_mat, varnames = c("Sig1", "Sig2"), value.name = "cor")

  # Preserve signature order
  melted$Sig1 <- factor(melted$Sig1, levels = sig_names)
  melted$Sig2 <- factor(melted$Sig2, levels = sig_names)

  ggplot(melted, aes(x = Sig2, y = Sig1, fill = cor)) +
    geom_tile(color = "grey30", linewidth = 0.3) +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0, limits = c(-1, 1),
      name = "Correlation",
      na.value = "white"
    ) +
    scale_y_discrete(limits = rev(sig_names)) +
    scale_x_discrete(position = "top") +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(
      axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 7),
      axis.text.y = element_text(size = 7),
      plot.title  = element_text(hjust = 0.5, size = 11),
      legend.position = "right"
    ) +
    coord_fixed()
}

# Output PDF in same directory as this script
out_pdf <- here::here(
  "code_for_internal_exploration/sig_correlation_heatmaps.pdf"
)
pdf(out_pdf, width = 12, height = 11)

# --- All samples ---
cors_all <- compute_cors(df)
print(make_tri_heatmap(cors_all$spearman,
                       "All samples - Spearman on log(counts + 1)"))
print(make_tri_heatmap(cors_all$pearson,
                       "All samples - Pearson on proportions"))

# --- Per cancer type ---
sample_ids <- colnames(df)
sample_cancer <- cancer_map[sample_ids]
cancer_types <- sort(unique(sample_cancer[!is.na(sample_cancer)]))

for (ct in cancer_types) {
  cols <- sample_ids[which(sample_cancer == ct)]
  n <- length(cols)
  if (n < 3) next  # need at least 3 samples for correlation

  sub <- df[, cols, drop = FALSE]
  cors_ct <- compute_cors(sub)

  print(make_tri_heatmap(
    cors_ct$spearman,
    paste0(ct, " (n=", n, ") - Spearman on log(counts + 1)")
  ))
  print(make_tri_heatmap(
    cors_ct$pearson,
    paste0(ct, " (n=", n, ") - Pearson on proportions")
  ))
}

dev.off()
message("Saved to ", out_pdf)
