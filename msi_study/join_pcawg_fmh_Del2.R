library(tidyverse)
library(openxlsx)
library(glue)

max_files = 200

pcawg <- read.csv(here::here(
  glue::glue("msi_study/check_pcawg_top{max_files}_Del2_U1_R5_9.csv")
))
fmh <- read.csv(here::here(
  glue::glue("msi_study/check_fmh_top{max_files}_Del2_U1_R5_9.csv")
))

# Extract total indel counts from ANY row, then remove it
pcawg_total <- pcawg$n[pcawg$short_visual == "ANY"]
fmh_total <- fmh$n[fmh$short_visual == "ANY"]

pcawg <- pcawg |>
  dplyr::filter(short_visual != "ANY") |>
  rename(pcawg_n = n)

fmh <- fmh |>
  dplyr::filter(short_visual != "ANY") |>
  rename(fmh_n = n)

joined <- full_join(pcawg, fmh, by = c("short_visual", "R")) |>
  mutate(
    pcawg_prop = pcawg_n / pcawg_total,
    fmh_prop = fmh_n / fmh_total,
    prop_ratio_pcawg_over_fmh = pcawg_prop / fmh_prop
  ) |>
  arrange(desc(R))

# Add ANY totals as the first row
any_row <- tibble(
  short_visual = "ANY",
  R = NA_integer_,
  pcawg_n = pcawg_total,
  fmh_n = fmh_total,
  pcawg_prop = NA_real_,
  fmh_prop = NA_real_,
  prop_ratio_pcawg_over_fmh = NA_real_
)
joined <- bind_rows(any_row, joined)

out_path <- here::here(glue::glue(
  "msi_study/joined_pcawg_fmh_top{max_files}_Del2_U1_R5_9.xlsx"
))
write.xlsx(joined, out_path)
cat(sprintf("Wrote %d rows to %s\n", nrow(joined), out_path))

# Join grouped by R only (summing across all short_visual)
pcawg_by_R <- pcawg |>
  group_by(R) |>
  summarise(pcawg_n = sum(pcawg_n), .groups = "drop")

fmh_by_R <- fmh |>
  group_by(R) |>
  summarise(fmh_n = sum(fmh_n), .groups = "drop")

joined_by_R <- full_join(pcawg_by_R, fmh_by_R, by = "R") |>
  mutate(
    pcawg_prop = pcawg_n / pcawg_total,
    fmh_prop = fmh_n / fmh_total,
    prop_ratio_pcawg_over_fmh = pcawg_prop / fmh_prop
  ) |>
  arrange(desc(R))

any_row_R <- tibble(
  R = NA_integer_,
  pcawg_n = pcawg_total,
  fmh_n = fmh_total,
  pcawg_prop = NA_real_,
  fmh_prop = NA_real_,
  prop_ratio_pcawg_over_fmh = NA_real_
)
joined_by_R <- bind_rows(any_row_R, joined_by_R)

out_path_R <- here::here(
  glue::glue(
    "msi_study/joined_pcawg_fmh_top{max_files}0_Del2_U1_R5_9_by_R.xlsx"
  )
)
write.xlsx(joined_by_R, out_path_R)
cat(sprintf("Wrote %d rows to %s\n", nrow(joined_by_R), out_path_R))
