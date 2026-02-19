library(tidyverse)
library(data.table)

# Patterns to track (same as msi_study.qmd tract-lengths block)
patterns <- c(
  "Del\\(T" = "Del(T)",
  "Ins\\(T" = "Ins(T)",
  "Del\\(C" = "Del(C)",
  "Ins\\(C" = "Ins(C)"
)

run_polyT_analysis <- function(vcf_dir, dataset_label, max_files, pdf_path) {
  all_files <- list.files(
    path.expand(vcf_dir),
    pattern = "annotated\\.indel\\.vcf\\.gz$",
    full.names = TRUE
  )
  file_sizes <- file.size(all_files)
  top_files <- all_files[order(file_sizes, decreasing = TRUE)][
    1:min(max_files, length(all_files))
  ]
  cat(sprintf("\n=== %s: selected %d files ===\n", dataset_label, length(top_files)))

  counts_list <- setNames(
    lapply(seq_along(patterns), function(x) tibble(repeat_count = integer(), n = integer())),
    names(patterns)
  )

  for (f in top_files) {
    cat("Reading", basename(f), "\n")
    vcf <- fread(f)

    for (i in seq_along(patterns)) {
      pat <- names(patterns)[i]

      counts <- vcf |>
        as_tibble() |>
        dplyr::filter(grepl(pat, Koh_476)) |>
        dplyr::count(R) |>
        dplyr::mutate(repeat_count = R)

      if (nrow(counts) > 0) {
        counts_list[[pat]] <- bind_rows(counts_list[[pat]], counts)
      }
    }
  }

  agg_list <- lapply(counts_list, function(tbl) {
    tbl |>
      group_by(repeat_count) |>
      summarise(n = sum(n), .groups = "drop")
  })

  pdf(pdf_path, width = 8, height = 5)
  for (i in seq_along(patterns)) {
    pat <- names(patterns)[i]
    label <- patterns[i]
    counts <- agg_list[[pat]]

    if (nrow(counts) == 0) {
      cat(sprintf("No mutations found for %s\n", label))
      next
    }

    cat(sprintf(
      "\n%s: Mean=%.2f, Median=%.1f, SD=%.2f, count=%d\n",
      label,
      weighted.mean(counts$repeat_count, counts$n),
      median(rep(counts$repeat_count, counts$n)),
      sd(rep(counts$repeat_count, counts$n)),
      sum(counts$n)
    ))

    p <- ggplot(counts, aes(x = repeat_count, weight = n)) +
      geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
      labs(
        title = paste(label, "repeat count distribution -",
                      length(top_files), dataset_label, "samples"),
        x = "Repeat count",
        y = "Count"
      ) +
      theme_minimal()
    print(p)
  }
  dev.off()
  cat(sprintf("\nPlots saved to %s\n", pdf_path))
}

# Run for PCAWG
run_polyT_analysis(
  vcf_dir = "~/MEGA/important_mut_sig_data/pcawg_indel_vcfs",
  dataset_label = "PCAWG",
  max_files = 1000,
  pdf_path = here::here("msi_study/checkSP_pcawg.pdf")
)

# Run for FMH
run_polyT_analysis(
  vcf_dir = "~/MEGA/important_mut_sig_data/fmh-unfiltered_vcfs",
  dataset_label = "FMH",
  max_files = 1000,
  pdf_path = here::here("msi_study/checkSP_fmh.pdf")
)
