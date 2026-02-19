library(tidyverse)
library(data.table)

# Patterns to track (same as msi_study.qmd tract-lengths block)
patterns <- c(
  "Del\\(T" = "Del(T)",
  "Ins\\(T" = "Ins(T)",
  "Del\\(C" = "Del(C)",
  "Ins\\(C" = "Ins(C)",
  "^Del2:U1:R\\(5,9\\)$" = "Del2:U1:R(5,9)"
)

max_files = 200

run_polyT_analysis <- function(
  vcf_dir,
  dataset_label,
  max_files,
  pdf_path,
  csv_path
) {
  all_files <- list.files(
    path.expand(vcf_dir),
    pattern = "annotated\\.indel\\.vcf\\.gz$",
    full.names = TRUE
  )
  file_sizes <- file.size(all_files)
  top_files <- all_files[order(file_sizes, decreasing = TRUE)][
    1:min(max_files, length(all_files))
  ]
  cat(sprintf(
    "\n=== %s: selected %d files ===\n",
    dataset_label,
    length(top_files)
  ))

  counts_list <- setNames(
    lapply(seq_along(patterns), function(x) {
      tibble(repeat_count = integer(), n = integer())
    }),
    names(patterns)
  )

  del2_detail <- tibble(
    short_visual = character(),
    R = integer(),
    n = integer()
  )
  del2_pat <- "^Del2:U1:R\\(5,9\\)$"
  total_indels <- 0L

  for (f in top_files) {
    cat("Reading", basename(f), "\n")
    vcf <- fread(f)
    total_indels <- total_indels + nrow(vcf)

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

    del2_rows <- vcf |>
      as_tibble() |>
      dplyr::filter(grepl(del2_pat, Koh_476)) |>
      dplyr::count(short_visual, R)
    if (nrow(del2_rows) > 0) {
      del2_detail <- bind_rows(del2_detail, del2_rows)
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

    n_lt9 <- sum(counts$n[counts$repeat_count < 9])
    n_ge9 <- sum(counts$n[counts$repeat_count >= 9])
    n_eq9 <- sum(counts$n[counts$repeat_count == 9])
    ratio_eq9_ge9 <- if (n_ge9 > 0) n_eq9 / n_ge9 else Inf

    cat(sprintf(
      "\n%s: Mean=%.2f, Median=%.1f, SD=%.2f, count=%d\n",
      label,
      weighted.mean(counts$repeat_count, counts$n),
      median(rep(counts$repeat_count, counts$n)),
      sd(rep(counts$repeat_count, counts$n)),
      sum(counts$n)
    ))
    cat(sprintf(
      "  n(R<9)=%d, n(R>=9)=%d, n(R==9)=%d, ratio(R==9 / R>=9)=%.4f\n",
      n_lt9, n_ge9, n_eq9, ratio_eq9_ge9
    ))

    p <- ggplot(counts, aes(x = repeat_count, weight = n)) +
      geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
      labs(
        title = paste(
          label,
          "repeat count distribution -",
          length(top_files),
          dataset_label,
          "samples"
        ),
        x = "Repeat count",
        y = "Count"
      ) +
      theme_minimal()
    print(p)
  }
  dev.off()
  cat(sprintf("\nPlots saved to %s\n", pdf_path))

  del2_agg <- del2_detail |>
    group_by(short_visual, R) |>
    summarise(n = sum(n), .groups = "drop") |>
    arrange(desc(n))
  any_row <- tibble(short_visual = "ANY", R = NA_integer_, n = total_indels)
  del2_agg <- bind_rows(any_row, del2_agg)
  write.csv(del2_agg, csv_path, row.names = FALSE)
  cat(sprintf(
    "Del2:U1:R(5,9) detail saved to %s (%d rows)\n",
    csv_path,
    nrow(del2_agg)
  ))
}

# Run for PCAWG - top 20
run_polyT_analysis(
  vcf_dir = "~/MEGA/important_mut_sig_data/pcawg_indel_vcfs",
  dataset_label = "PCAWG",
  max_files = max_files,
  pdf_path = here::here(glue::glue("msi_study/check_pcawg_top{max_files}.pdf")),
  csv_path = here::here(glue::glue(
    "msi_study/check_pcawg_top{max_files}_Del2_U1_R5_9.csv"
  ))
)

# Run for FMH - top 20
run_polyT_analysis(
  vcf_dir = "~/MEGA/important_mut_sig_data/fmh-unfiltered_vcfs",
  dataset_label = "FMH",
  max_files = max_files,
  pdf_path = here::here(glue::glue("msi_study/check_fmh_top{max_files}.pdf")),
  csv_path = here::here(glue::glue(
    "msi_study/check_fmh_top{max_files}_Del2_U1_R5_9.csv"
  ))
)
