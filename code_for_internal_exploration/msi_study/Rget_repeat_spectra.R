library(tidyverse)
library(data.table)
library(openxlsx)

get_repeat_spectra <- function(
  sample_id_list,
  dataset_label,
  pdf_path,
  output_csv_path
) {
  # Patterns to track (same as msi_study.qmd tract-lengths block)
  patterns <- c(
    "Del\\(T" = "Del(T)",
    "Ins\\(T" = "Ins(T)",
    "Del\\(C" = "Del(C)",
    "Ins\\(C" = "Ins(C)",
    "^Del2:U1:R\\(5,9\\)$" = "Del2:U1:R(5,9)"
  )

  if (grepl("::SP", sample_id_list[1], fixed = TRUE)) {
    vcf_dir <- "~/MEGA/important_mut_sig_data/pcawg_indel_vcfs"
  } else {
    vcf_dir <- "~/MEGA/important_mut_sig_data/fmh-unfiltered_vcfs"
  }

  all_files <- list.files(
    path.expand(vcf_dir),
    pattern = "annotated\\.indel\\.vcf\\.gz$",
    full.names = TRUE
  )

  # Extract sample IDs from filenames (strip tissue prefix)
  sample_ids_from_ids <- sub("^.*::", "", sample_id_list)
  # Match files by sample ID in basename
  matched_files <- all_files[
    sub("\\.annotated\\.indel\\.vcf\\.gz$", "", basename(all_files)) %in%
      sample_ids_from_ids
  ]

  cat(sprintf(
    "\n=== %s: matched %d of %d requested samples ===\n",
    dataset_label,
    length(matched_files),
    length(sample_id_list)
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

  for (f in matched_files) {
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

  summary_dt <- data.table(
    pattern = character(),
    mean = numeric(),
    median = numeric(),
    sd = numeric(),
    count = integer(),
    n_R_lt9 = integer(),
    n_R_ge9 = integer(),
    n_R_eq9 = integer(),
    ratio_eq9_ge9 = numeric()
  )

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

    summary_dt <- rbind(
      summary_dt,
      data.table(
        pattern = label,
        mean = weighted.mean(counts$repeat_count, counts$n),
        median = median(rep(counts$repeat_count, counts$n)),
        sd = sd(rep(counts$repeat_count, counts$n)),
        count = sum(counts$n),
        n_R_lt9 = n_lt9,
        n_R_ge9 = n_ge9,
        n_R_eq9 = n_eq9,
        ratio_eq9_ge9 = ratio_eq9_ge9
      )
    )

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
      n_lt9,
      n_ge9,
      n_eq9,
      ratio_eq9_ge9
    ))

    p <- ggplot(counts, aes(x = repeat_count, weight = n)) +
      geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
      labs(
        title = paste(
          label,
          "repeat count distribution -",
          length(matched_files),
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

  xlsx_path <- here::here(glue::glue(
    "msi_study/{dataset_label}_proportion_summary.xlsx"
  ))
  write.xlsx(summary_dt, xlsx_path)
  cat(sprintf("Proportion summary saved to %s\n", xlsx_path))

  del2_agg <- del2_detail |>
    group_by(short_visual, R) |>
    summarise(n = sum(n), .groups = "drop") |>
    arrange(desc(n))
  any_row <- tibble(short_visual = "ANY", R = NA_integer_, n = total_indels)
  del2_agg <- bind_rows(any_row, del2_agg)
  write.csv(del2_agg, output_csv_path, row.names = FALSE)
  cat(sprintf(
    "Del2:U1:R(5,9) detail saved to %s (%d rows)\n",
    output_csv_path,
    nrow(del2_agg)
  ))
}

# Read sample lists and fix ".." -> "::"
hits_7 <- read.csv(here::here("msi_study/7_hits_ge_0.0.csv"))
hits_7$spectrum <- gsub("..", "::", hits_7$spectrum, fixed = TRUE)

hits_J <- read.csv(here::here("msi_study/J_hits_ge_0.0.csv"))
hits_J$spectrum <- gsub("..", "::", hits_J$spectrum, fixed = TRUE)

# Run for PCAWG_7
get_repeat_spectra(
  sample_id_list = hits_7$spectrum,
  dataset_label = "PCAWG_7",
  pdf_path = here::here("msi_study/check_PCAWG_7.pdf"),
  output_csv_path = here::here("msi_study/check_PCAWG_7_Del2_U1_R5_9.csv")
)

# Run for FMH_J
get_repeat_spectra(
  sample_id_list = hits_J$spectrum,
  dataset_label = "FMH_J",
  pdf_path = here::here("msi_study/check_FMH_J.pdf"),
  output_csv_path = here::here("msi_study/check_FMH_J_Del2_U1_R5_9.csv")
)
