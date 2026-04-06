# Plot histograms of L, mh_length, and mutations per tumor from deep_dive_H.csv
# (R_intuitive is in separate scripts)
#
# Skin vs non-Skin comparison: each page has Skin on top, non-Skin below,
# with shared x-axis limits and ticks.
#   Page 1: L (indel length)
#   Page 2: Mutations per tumor
#   Page 3: mh_length
#
# Output:
#   code_for_internal_exploration/indel_deep_dive/deep_dive_H_histograms.pdf
#
# Usage:
#   Rscript plot_deep_dive_H.R

library(ggplot2)
library(gridExtra)
library(scales)
library(here)

integer_breaks <- function(x) {
  lo <- floor(min(x))
  hi <- ceiling(max(x))
  rng <- hi - lo
  by <- max(1, round(rng / 6))
  seq(lo, hi, by = by)
}

# Like integer_breaks but always starts at 0 (for y-axes / counts)
y_breaks <- function(x) {
  hi <- ceiling(max(x))
  by <- max(1, round(hi / 6))
  seq(0, hi, by = by)
}

df <- read.csv(
  here::here("code_for_internal_exploration/indel_deep_dive/deep_dive_H.csv"),
  check.names = FALSE
)

out_path <- here::here(
  "code_for_internal_exploration/indel_deep_dive",
  "deep_dive_H_histograms.pdf"
)

min_mutations <- 11

# Filter to tumors with >= min_mutations
muts_all <- aggregate(count ~ tumor_id, data = df, FUN = sum)
keep_tumors <- muts_all$tumor_id[muts_all$count >= min_mutations]
df <- df[df$tumor_id %in% keep_tumors, ]

# Split into Skin and non-Skin
skin <- df[df$cancer_type == "Skin", ]
non_skin <- df[df$cancer_type != "Skin", ]

# Expand by count
skin_exp <- skin[rep(seq_len(nrow(skin)), skin$count), ]
ns_exp <- non_skin[rep(seq_len(nrow(non_skin)), non_skin$count), ]

n_skin <- length(unique(skin$tumor_id))
n_ns <- length(unique(non_skin$tumor_id))
n_ns_types <- length(unique(non_skin$cancer_type))

sub_skin <- paste0("Skin: ", n_skin, " tumors")
sub_ns <- paste0(
  "Non-Skin: ", n_ns, " tumors across ", n_ns_types, " cancer types"
)

# Mutations per tumor
muts_skin <- muts_all[muts_all$tumor_id %in% unique(skin$tumor_id), ]
muts_ns <- muts_all[muts_all$tumor_id %in% unique(non_skin$tumor_id), ]
muts_ns <- merge(muts_ns, unique(non_skin[, c("tumor_id", "cancer_type")]),
                 by = "tumor_id")

pdf(out_path, width = 8.5, height = 11)

# --- Page 1: L ---
all_L <- c(skin_exp$L[!is.na(skin_exp$L)], ns_exp$L[!is.na(ns_exp$L)])
L_breaks <- integer_breaks(all_L)
L_lim <- range(L_breaks)

p_L_skin <- ggplot(skin_exp[!is.na(skin_exp$L), ], aes(x = L)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
  scale_x_continuous(breaks = L_breaks, limits = L_lim) +
  scale_y_continuous(breaks = y_breaks) +
  labs(title = "L (indel length) — Skin", subtitle = sub_skin,
       x = "L", y = "Number of mutations") +
  theme_minimal()

p_L_ns <- ggplot(ns_exp[!is.na(ns_exp$L), ], aes(x = L, fill = cancer_type)) +
  geom_histogram(binwidth = 1, color = "white") +
  scale_x_continuous(breaks = L_breaks, limits = L_lim) +
  scale_y_continuous(breaks = y_breaks) +
  labs(title = "L (indel length) — Non-Skin", subtitle = sub_ns,
       x = "L", y = "Number of mutations", fill = "Cancer type") +
  theme_minimal()

grid.arrange(p_L_skin, p_L_ns, ncol = 1)

# --- Page 2: Mutations per tumor ---
all_muts <- c(muts_skin$count, muts_ns$count)
M_breaks <- integer_breaks(all_muts)
M_lim <- range(M_breaks)
M_bw <- max(1, round(diff(range(all_muts)) / 20))

p_M_skin <- ggplot(muts_skin, aes(x = count)) +
  geom_histogram(binwidth = M_bw, fill = "purple", color = "white") +
  scale_x_continuous(breaks = M_breaks, limits = M_lim) +
  scale_y_continuous(breaks = y_breaks) +
  labs(title = "Mutations per tumor — Skin", subtitle = sub_skin,
       x = "Mutations per tumor", y = "Number of tumors") +
  theme_minimal()

p_M_ns <- ggplot(muts_ns, aes(x = count, fill = cancer_type)) +
  geom_histogram(binwidth = M_bw, color = "white") +
  scale_x_continuous(breaks = M_breaks, limits = M_lim) +
  scale_y_continuous(breaks = y_breaks) +
  labs(title = "Mutations per tumor — Non-Skin", subtitle = sub_ns,
       x = "Mutations per tumor", y = "Number of tumors", fill = "Cancer type") +
  theme_minimal()

grid.arrange(p_M_skin, p_M_ns, ncol = 1)

# --- Page 3: mh_length ---
all_mh <- c(
  skin_exp$mh_length[!is.na(skin_exp$mh_length)],
  ns_exp$mh_length[!is.na(ns_exp$mh_length)]
)
mh_breaks <- integer_breaks(all_mh)
mh_lim <- range(mh_breaks)

p_mh_skin <- ggplot(
  skin_exp[!is.na(skin_exp$mh_length), ], aes(x = mh_length)
) +
  geom_histogram(binwidth = 1, fill = "darkgreen", color = "white") +
  scale_x_continuous(breaks = mh_breaks, limits = mh_lim) +
  scale_y_continuous(breaks = y_breaks) +
  labs(title = "mh_length — Skin", subtitle = sub_skin,
       x = "mh_length", y = "Number of mutations") +
  theme_minimal()

p_mh_ns <- ggplot(
  ns_exp[!is.na(ns_exp$mh_length), ],
  aes(x = mh_length, fill = cancer_type)
) +
  geom_histogram(binwidth = 1, color = "white") +
  scale_x_continuous(breaks = mh_breaks, limits = mh_lim) +
  scale_y_continuous(breaks = y_breaks) +
  labs(title = "mh_length — Non-Skin", subtitle = sub_ns,
       x = "mh_length", y = "Number of mutations", fill = "Cancer type") +
  theme_minimal()

grid.arrange(p_mh_skin, p_mh_ns, ncol = 1)

dev.off()
message("Wrote ", out_path)
