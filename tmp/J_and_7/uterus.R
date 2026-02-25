library(data.table)
library(dplyr)
library(tidyverse)
to.test = data.table::fread(here::here(
  "tmp/J_and_7/DRUP01010022T.annotated.indels.txt"
))

to.test %>% filter(is.na(Koh_code_89)) -> to.test.na
unique(to.test.na$Koh_89)
unique(to.test.na$Koh_476)

# calculate it a different way
to.test %>% filter(Koh_476 == "Del2:U1:R(5,9)") %>% filter(is.na(Koh_code_89))

to.test %>% filter(Koh_476 == "Del2:U1:R(5,9)") -> mmdel2

newannot = fread(here::here("tmp/J_and_7/DRUP01010022T.annotated.indel.vcf.gz"))
newannot %>% filter(Koh_476 == "Del2:U1:R(5,9)") -> newdel2
table(newdel2$R)
