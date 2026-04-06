# R_intuitive stacked histogram for non-Skin tumors with InsDel_H activity

library(ggplot2)
library(here)

df <- read.csv(
  here::here("code_for_internal_exploration/indel_deep_dive/deep_dive_H.csv")
)

ns <- df[df$cancer_type != "Skin" & !is.na(df$R_intuitive), ]
ns_exp <- ns[rep(seq_len(nrow(ns)), ns$count), ]
n_tumors <- length(unique(ns$tumor_id))
n_types <- length(unique(ns$cancer_type))

p <- ggplot(ns_exp, aes(x = R_intuitive, fill = cancer_type)) +
  geom_histogram(binwidth = 1, color = "white") +
  labs(
    title = "R_intuitive — Non-Skin",
    subtitle = paste0(n_tumors, " tumors across ", n_types, " cancer types, ",
                      nrow(ns_exp), " mutations"),
    x = "R_intuitive",
    y = "Number of mutations",
    fill = "Cancer type"
  ) +
  theme_minimal()

out <- here::here(
  "code_for_internal_exploration/indel_deep_dive",
  "R_intuitive_non_skin.pdf"
)
pdf(out, width = 8, height = 5)
print(p)
dev.off()
message("Wrote ", out)
