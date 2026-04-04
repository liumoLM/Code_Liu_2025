#!/usr/bin/env Rscript

# Read the 89-type signature assignment matrix (signatures x samples),
# transpose it, join with cancer type, and produce a summary table:
# cancer_type | number_of_samples | <sig1> | <sig2> | ...
# where each signature column counts how many samples have that signature (value > 0).

library(tidyverse)

assignment <- read.delim(
  here::here("Manuscript_data/finalized_cap9/liu_et_al_89_assignment.tsv"),
  row.names   = 1,
  check.names = FALSE
)

# Transpose: samples become rows, signatures become columns
assignment_t <- as.data.frame(t(assignment))
assignment_t$sampleid <- rownames(assignment_t)

sample_info <- read.delim(
  here::here("Manuscript_data/sample_info.tsv")
)

# Join on Patient = sampleid
merged <- inner_join(
  sample_info[, c("Patient", "Cancer_Type")],
  assignment_t,
  by = c("Patient" = "sampleid")
)

sig_names <- colnames(assignment_t)[colnames(assignment_t) != "sampleid"]

summary_table <- merged %>%
  group_by(Cancer_Type) %>%
  summarise(
    number_of_samples = n(),
    across(
      all_of(sig_names),
      ~ sum(. > 0)
    ),
    .groups = "drop"
  ) %>%
  rename(cancer_type = Cancer_Type)

# Rename signature columns to just the signature name (they already are)
write_tsv(summary_table, here::here("plot_output/samples_with_sig_by_cancer_type.tsv"))
cat("Wrote", nrow(summary_table), "cancer types x", length(sig_names), "signatures\n")
cat("Output: plot_output/samples_with_sig_by_cancer_type.tsv\n")
