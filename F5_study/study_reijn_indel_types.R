# Study specifiy mutations types in the Reijn et al data

library(magrittr)

# source("code/generate_89_and_476_type_catalog.R")
source("code/Generate_Koh89_Koh476_catalog_0121.R")

files = dir("Manuscript_data", pattern = "Reijn", full.names = TRUE)
mouse = files[1]
rpe1 = files[2]

mvcf = read.delim(mouse, sep = '\t')

rvcf = read.delim(rpe1, sep = '\t')

mvcf %>%
  dplyr::filter(nchar(ins_or_del_seq) > 1) %>%
  dplyr::count(short_visual) %>%
  dplyr::arrange(desc(n)) %>%
  write.csv("F5_study/reijn_mouse.csv")

rvcf %>%
  dplyr::filter(nchar(ins_or_del_seq) > 1) %>%
  dplyr::count(short_visual) %>%
  dplyr::arrange(desc(n)) %>%
  write.csv("F5_study/reijn_rpe1.csv")
