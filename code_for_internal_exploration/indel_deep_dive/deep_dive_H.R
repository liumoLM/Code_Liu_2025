# Deep dive into InsDel_H: examine all tumors with InsDel_H activity
#
# Looks for Del(6,):U(4,):R(2,9) and Del(7,):M(6,) patterns in annotated VCFs
# of all tumors with any InsDel_H mutations.
#
# Usage:
#   Rscript deep_dive_H.R

library(here)

source(here::here("code_for_internal_exploration/indel_deep_dive.R"))

# Load assignment matrix and sample info
assignments <- read.delim(
  here::here("Manuscript_data/finalized_cap9/liu_et_al_89_assignment.tsv"),
  row.names = 1,
  check.names = FALSE
)

sample_info <- read.delim(
  here::here("Manuscript_data/sample_info.tsv"),
  check.names = FALSE
)

# Find all tumors with any mutations attributed to InsDel_H
h_values <- as.numeric(assignments["InsDel_H", ])
names(h_values) <- colnames(assignments)
h_values <- sort(h_values[h_values > 0], decreasing = TRUE)

cat("Tumors with InsDel_H activity:", length(h_values), "\n")
cat("Top 10:\n")
print(head(h_values, 10))

# Patterns of interest for InsDel_H
# Del(6,):U(4,):R(2,9) — long deletions with 4+ unit repeats
# Del(7,):M(6,)        — long deletions with 6+ microhomology bases
patterns_to_match <- c(
  "^Del\\(6,\\):U\\(4,\\):R\\(2,9\\)$" = "Del(6,):U(4,):R(2,9)",
  "^Del\\(7,\\):M\\(6,\\)$" = "Del(7,):M(6,)"
)

result <- indel_deep_dive(
  samples_to_fetch = names(h_values),
  patterns_to_match = patterns_to_match,
  sample_info = sample_info,
  cap9 = FALSE
)

csv_path <- here::here(
  "code_for_internal_exploration/indel_deep_dive",
  paste0("deep_dive_H.csv")
)
write.csv(result, csv_path, row.names = FALSE)
cat("Saved visual table to", csv_path, "\n")
