# R_intuitive histogram for Skin tumors with InsDel_H activity

library(ggplot2)
library(here)

df <- read.csv(
  here::here("code_for_internal_exploration/indel_deep_dive/deep_dive_H.csv")
)

skin <- df[df$cancer_type == "Skin" & !is.na(df$R_intuitive), ]
skin_exp <- skin[rep(seq_len(nrow(skin)), skin$count), ]
n_tumors <- length(unique(skin$tumor_id))

p <- ggplot(skin_exp, aes(x = R_intuitive)) +
  geom_histogram(binwidth = 1, fill = "darkorange", color = "white") +
  labs(
    title = "R_intuitive — Skin",
    subtitle = paste0(n_tumors, " tumors, ", nrow(skin_exp), " mutations"),
    x = "R_intuitive",
    y = "Number of mutations"
  ) +
  theme_minimal()

out <- here::here(
  "code_for_internal_exploration/indel_deep_dive",
  "R_intuitive_skin.pdf"
)
pdf(out, width = 8, height = 5)
print(p)
dev.off()
message("Wrote ", out)
