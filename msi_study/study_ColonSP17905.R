library(tidyr)
library(dplyr)

sample_id = "Colon::SP17905"

source(here::here("code/read_annotated_vcf.R"))

vv = read_annotated_vcf(sample_id)

vv %>%
  count(Koh_476, R, sort = TRUE, name = "newcount") -> vvv


c476 = ICAMS::annot_vcf_to_476_catalog(vv)
c89 = ICAMS::annot_vcf_to_89_catalog(vv)
c83 = ICAMS::annot_vcf_to_83_catalog(vv)


id = "Pancreas::SP125746"
vv2 = read_annotated_vcf(id)
