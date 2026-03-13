#!/usr/bin/env Rscript
# Render cluster_cap9_koh89.qmd with a parameterized min_similarity
# and open the resulting HTML file.
#
# Usage:
#   Rscript render_cluster.R --min-similarity 0.95

library(argparser)

p <- arg_parser("Render cluster_cap9_koh89.qmd with parameterized min_similarity")
p <- add_argument(p, "--min-similarity", type = "numeric", default = 0.95,
                  help = "Minimum cosine similarity for clustering")
args <- parse_args(p)

script_dir  <- here::here("Manuscript_data", "Mo_CAP9_analysis")
qmd_file    <- file.path(script_dir, "cluster_cap9_koh89.qmd")
output_file <- sprintf("cluster_cap9_89_minsim_%s.html", args$min_similarity)

message("Rendering with min_similarity = ", args$min_similarity)
message("Output: ", output_file)

quarto::quarto_render(
  input         = qmd_file,
  output_file   = output_file,
  execute_params = list(min_similarity = args$min_similarity)
)

output_path <- file.path(script_dir, output_file)
message("Opening ", output_path)
system2("xdg-open", output_path, wait = FALSE,
        stdout = FALSE, stderr = FALSE)
