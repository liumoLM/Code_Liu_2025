library(tidyr)
library(dplyr)

sample_id = "Colon::SP17905"

source(here::here("code/read_annotated_vcf.R"))

vv = read_annotated_vcf(sample_id)

vv %>%
  count(Koh_476, R, sort = TRUE, name = "newcount") -> vvv
