library(tidyverse)

pcawg <- read.csv(here::here("msi_study/check_pcawg_Del2_U1_R5_9.csv")) |>
  rename(pcawg_R = R, pcawg_n = n)

fmh <- read.csv(here::here("msi_study/check_fmh_Del2_U1_R5_9.csv")) |>
  rename(fmh_R = R, fmh_n = n)

joined <- full_join(pcawg, fmh, by = "short_visual",
                    relationship = "many-to-many") |>
  arrange(desc(fmh_n))

out_path <- here::here("msi_study/joined_pcawg_fmh_Del2_U1_R5_9.csv")
write.csv(joined, out_path, row.names = FALSE)
cat(sprintf("Wrote %d rows to %s\n", nrow(joined), out_path))
