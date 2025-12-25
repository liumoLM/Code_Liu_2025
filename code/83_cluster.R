library(ggplot2)
library(ggdendro)

# Read the data
cosmic_sigs <- read.table(
  "data/COSMIC_v3.5_ID_GRCh37.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1
)

our_sigs <- read.table(
  "data/type83_our_sigs.tsv",
  header = TRUE,
  sep = "\t",
  row.names = 1
)

jin_sigs <- read.table(
  "data/jin_2024_indel_sigs_sup_tab_1.tsv",
  header = TRUE,
  sep = "\t",
  row.names = 1
)

# Combine and transpose (samples as rows, features as columns)
combined <- cbind(cosmic_sigs, our_sigs, jin_sigs)
combined_t <- t(combined)

# Normalize rows to unit vectors for cosine distance
combined_norm <- combined_t / sqrt(rowSums(combined_t^2))

# Compute cosine distance matrix (1 - cosine similarity)
cosine_sim <- combined_norm %*% t(combined_norm)
cosine_dist <- as.dist(1 - cosine_sim)

# Hierarchical clustering (average linkage works well with cosine distance)
hc <- hclust(cosine_dist, method = "average")

# Create source vector for coloring
source_vec <- c(
  rep("COSMIC", ncol(cosmic_sigs)),
  rep("Ours", ncol(our_sigs)),
  rep("Jin", ncol(jin_sigs))
)
names(source_vec) <- rownames(combined_t)

# Extract dendrogram data
dend_data <- dendro_data(hc)

# Add source info to labels
label_df <- dend_data$labels
label_df$source <- source_vec[label_df$label]

# Color-blind friendly palette (Okabe-Ito)
cb_palette <- c(
  "COSMIC" = "#0072B2",
  "Ours" = "#D55E00",
  "Jin" = "#009E73"
)

# Plot
p <- ggplot() +
  geom_segment(
    data = dend_data$segments,
    aes(x = x, y = y, xend = xend, yend = yend)
  ) +
  geom_text(
    data = label_df,
    aes(x = x, y = y - 0.01, label = label, color = source),
    hjust = 1,
    angle = 90,
    size = 3
  ) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = cb_palette) +
  coord_cartesian(ylim = c(-0.15, NA), clip = "off") +
  labs(
    title = "Hierarchical Clustering of ID83 Signatures (Cosine Distance)",
    x = "",
    y = "Distance",
    color = "Source"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 60, l = 10)
  )

# Save to PDF
cairo_pdf("hclust_id83_signatures.pdf", width = 16, height = 8)
print(p)
dev.off()

message("Created hclust_id83_signatures.pdf")
