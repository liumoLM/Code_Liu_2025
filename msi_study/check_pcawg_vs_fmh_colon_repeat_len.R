library(tidyverse)
library(data.table)
library(here)
library(openxlsx)

# Patterns to track (same as msi_study.qmd tract-lengths block)
patterns <- c(
  "Del\\(T" = "Del(T)",
  "Ins\\(T" = "Ins(T)",
  "Del\\(C" = "Del(C)",
  "Ins\\(C" = "Ins(C)",
  "^Del2:U1:R\\(5,9\\)$" = "Del2:U1:R(5,9)"
)

mincolsum = 14000

run_analysis <- function(
  vcf_dir,
  dataset_label,
  pdf_path,
  csv_path
) {
  colonsamples = data.table::fread(here(
    "Manuscript_data/cancertype_to_sampleid.csv"
  )) %>%
    dplyr::filter(cancertype %in% c("Colon", "Uterus", "Prostate")) %>%
    pull(sampleid)

  # Filter colonsamples to those with total indel count >= mincolsum
  spectra <- read.delim(
    here("Manuscript_data/Liu_et_al_83_type_spectra.tsv"),
    row.names = 1,
    check.names = FALSE
  )
  col_sampleids <- sub("^.*::", "", colnames(spectra))
  colon_idx <- which(col_sampleids %in% colonsamples)
  colon_colsums <- colSums(spectra[, colon_idx])
  colonsamples <- col_sampleids[colon_idx][colon_colsums >= mincolsum]
  cat(sprintf(
    "%s: %d colon samples with colSums >= %d\n",
    dataset_label,
    length(colonsamples),
    mincolsum
  ))
  # browser()
  candidate_files <- file.path(
    path.expand(vcf_dir),
    paste0(colonsamples, ".annotated.indel.vcf.gz")
  )
  top_files <- candidate_files[file.exists(candidate_files)]
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

  summary_rows <- list()

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

    summary_rows[[length(summary_rows) + 1]] <- tibble(
      mutation_type = label,
      mean_R = weighted.mean(counts$repeat_count, counts$n),
      median_R = median(rep(counts$repeat_count, counts$n)),
      sd_R = sd(rep(counts$repeat_count, counts$n)),
      count_indels = sum(counts$n),
      `n(R<=9)` = n_lt9 + n_eq9,
      `n(R>9)` = n_ge9 - n_eq9,
      `n(R<9)` = n_lt9,
      `n(R>=9)` = n_ge9,
      `n(R==9)` = n_eq9,
      `ratio(R==9/R>=9)` = ratio_eq9_ge9
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

  summary_df <- bind_rows(summary_rows)

  # Round for display: 3 decimals for mean_R, sd_R, ratio(R==9/R>=9)
  summary_df$mean_R <- round(summary_df$mean_R, 3)
  summary_df$sd_R <- round(summary_df$sd_R, 3)
  summary_df$`ratio(R==9/R>=9)` <- round(summary_df$`ratio(R==9/R>=9)`, 3)

  # Write with openxlsx for bold header row
  wb <- createWorkbook()
  addWorksheet(wb, "Summary")

  # Row 1: dataset_label in bold
  bold_style <- createStyle(textDecoration = "bold")
  writeData(wb, "Summary", data.frame(x = dataset_label), startRow = 1,
            colNames = FALSE)
  addStyle(wb, "Summary", bold_style, rows = 1, cols = 1)

  # Row 2 onward: summary_df with header (no ratio rows yet)
  writeData(wb, "Summary", summary_df, startRow = 2)
  addStyle(wb, "Summary", bold_style, rows = 2, cols = seq_along(summary_df))

  # Excel rows: header=2, Del(T)=3, Ins(T)=4, Del(C)=5, Ins(C)=6, Del2=7
  # Empty row = 8, ratio1 = 9, ratio2 = 10
  # Columns: E=count_indels, F=n(R<=9), G=n(R>9)
  n_data <- nrow(summary_df)
  empty_row <- n_data + 3  # row after last data row
  r1 <- empty_row + 1
  r2 <- empty_row + 2
  del2_excel_row <- 2 + which(summary_df$mutation_type == "Del2:U1:R(5,9)")
  delT_excel_row <- 2 + which(summary_df$mutation_type == "Del(T)")
  mono_rows <- 2 + which(
    summary_df$mutation_type %in% c("Del(T)", "Ins(T)", "Del(C)", "Ins(C)")
  )

  # Ratio row labels
  writeData(wb, "Summary", "Del2:U1:R(5,9) / Del(T)",
            startRow = r1, startCol = 1, colNames = FALSE)
  writeData(wb, "Summary", "Del2:U1:R(5,9) / (Del(T)+Ins(T)+Del(C)+Ins(C))",
            startRow = r2, startCol = 1, colNames = FALSE)

  # 4-decimal format for ratio cells
  fmt4 <- createStyle(numFmt = "0.0000")

  # Ratio1: Del2 / Del(T), denominator for F,G is count_indels (col E)
  # E: count_indels / count_indels
  writeFormula(wb, "Summary", startRow = r1, startCol = 5,
    x = sprintf("=E%d/E%d", del2_excel_row, delT_excel_row))
  # F: n(R<=9) / Del(T) count_indels
  writeFormula(wb, "Summary", startRow = r1, startCol = 6,
    x = sprintf("=F%d/E%d", del2_excel_row, delT_excel_row))
  # G: n(R>9) / Del(T) count_indels
  writeFormula(wb, "Summary", startRow = r1, startCol = 7,
    x = sprintf("=G%d/E%d", del2_excel_row, delT_excel_row))

  addStyle(wb, "Summary", fmt4, rows = r1, cols = 5:7)

  # Ratio2: Del2 / sum(mono), denominator for F,G is sum of count_indels (col E)
  mono_sum_E <- paste0("E", mono_rows, collapse = "+")
  writeFormula(wb, "Summary", startRow = r2, startCol = 5,
    x = sprintf("=E%d/(%s)", del2_excel_row, mono_sum_E))
  writeFormula(wb, "Summary", startRow = r2, startCol = 6,
    x = sprintf("=F%d/(%s)", del2_excel_row, mono_sum_E))
  writeFormula(wb, "Summary", startRow = r2, startCol = 7,
    x = sprintf("=G%d/(%s)", del2_excel_row, mono_sum_E))

  addStyle(wb, "Summary", fmt4, rows = r2, cols = 5:7)

  # Set column widths to at least the width of the header
  col_widths <- pmax(nchar(colnames(summary_df)) + 2, 8)
  setColWidths(wb, "Summary", cols = seq_along(summary_df), widths = col_widths)

  xlsx_path <- here::here(
    "msi_study",
    paste0("check_", dataset_label, "_summary.xlsx")
  )
  saveWorkbook(wb, xlsx_path, overwrite = TRUE)
  cat(sprintf("Summary saved to %s\n", xlsx_path))

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

# Run for PCAWG
run_analysis(
  vcf_dir = "~/MEGA/important_mut_sig_data/pcawg_indel_vcfs",
  dataset_label = "PCAWG",
  pdf_path = here::here("msi_study/check_pcawg_colon_etc.pdf"),
  csv_path = here::here("msi_study/check_pcawg_colon_etc_Del2_U1_R5_9.csv")
)

# Run for FMH
run_analysis(
  vcf_dir = "~/MEGA/important_mut_sig_data/fmh-unfiltered_vcfs",
  dataset_label = "FMH",
  pdf_path = here::here("msi_study/check_fmh_etc_colon.pdf"),
  csv_path = here::here("msi_study/check_fmh_etc_Del2_U1_R5_9.csv")
)


run_analysis(
  vcf_dir = "~/MEGA/important_mut_sig_data/ICGC-Pan-Can-PCAWG7-2016-08-12-rawdata/final_consensus_12aug_passonly/graylist/indel",
  dataset_label = "PCAWG_gray",
  pdf_path = here::here("msi_study/check_pcawg_gray.pdf"),
  csv_path = here::here("msi_study/check_pcawg_gray_Del2_U1_R5_9.csv")
)
