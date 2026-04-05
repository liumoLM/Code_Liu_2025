library(tidyverse)
library(data.table)
library(here)

source(here::here("code_for_internal_exploration/read_annotated_vcf.R"))

#' Analyze indel repeat-count distributions from annotated VCFs
#'
#' For each pattern, aggregates repeat-count distributions across the
#' requested samples and produces histograms (PDF). Returns a table of
#' short_visual counts per pattern.
#'
#' @param samples_to_fetch Character vector of tumor IDs (optionally with
#'   "CancerType::" prefix). Each is fetched via \code{read_annotated_vcf()}.
#' @param patterns_to_match Named character vector: regex names (matched against
#'   Koh_476 column) mapped to display labels.
#' @param sample_info Data frame with columns Patient and Cancer_Type for
#'   annotating the output. If NULL, cancer_type column will be NA.
#' @param cap9 Logical: if TRUE, filter to R <= 9 before counting.
#' @param pdf_path Output PDF path. Defaults to a timestamped file in
#'   code_for_internal_exploration/indel_deep_dive/.
#' @return A tibble with columns: tumor_id, cancer_type, pattern,
#'   short_visual, R, L, count.
indel_deep_dive <- function(
  samples_to_fetch,
  patterns_to_match,
  sample_info = NULL,
  cap9 = TRUE,
  pdf_path = NULL
) {
  out_dir <- here::here("code_for_internal_exploration/indel_deep_dive")
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

  if (is.null(pdf_path)) {
    pdf_path <- file.path(
      out_dir,
      paste0("indel_deep_dive_", timestamp, ".pdf")
    )
  }

  counts_list <- setNames(
    lapply(seq_along(patterns_to_match), function(x) {
      tibble(repeat_count = integer(), n = integer())
    }),
    names(patterns_to_match)
  )

  # Collect short_visual counts per pattern
  visual_list <- list()

  total_indels <- 0L
  n_loaded <- 0L

  for (sample_id in samples_to_fetch) {
    cat("Fetching", sample_id, "\n")
    vcf <- read_annotated_vcf(sample_id)
    if (is.null(vcf)) {
      cat("  Skipped (not found)\n")
      next
    }
    n_loaded <- n_loaded + 1L
    total_indels <- total_indels + nrow(vcf)

    if (cap9) {
      vcf <- dplyr::filter(vcf, R <= 9)
    }

    for (i in seq_along(patterns_to_match)) {
      pat <- names(patterns_to_match)[i]

      counts <- vcf |>
        as_tibble() |>
        dplyr::filter(grepl(pat, Koh_476)) |>
        dplyr::count(R) |>
        dplyr::mutate(repeat_count = R)

      if (nrow(counts) > 0) {
        counts_list[[pat]] <- bind_rows(counts_list[[pat]], counts)
      }

      # Collect short_visual, R, L, and prefix counts for this pattern
      vis_counts <- vcf |>
        as_tibble() |>
        dplyr::filter(grepl(pat, Koh_476)) |>
        dplyr::mutate(prefix = sub(".*([A-Z]) <.*", "\\1", long_visual)) |>
        dplyr::count(short_visual, R, L, prefix)
      if (nrow(vis_counts) > 0) {
        vis_counts$pattern <- patterns_to_match[i]
        vis_counts$tumor_id <- sample_id
        visual_list <- c(visual_list, list(vis_counts))
      }
    }

  }

  cat(sprintf(
    "\nLoaded %d of %d requested samples\n",
    n_loaded,
    length(samples_to_fetch)
  ))

  agg_list <- lapply(counts_list, function(tbl) {
    tbl |>
      group_by(repeat_count) |>
      summarise(n = sum(n), .groups = "drop")
  })

  pdf(pdf_path, width = 8, height = 5)
  for (i in seq_along(patterns_to_match)) {
    pat <- names(patterns_to_match)[i]
    label <- patterns_to_match[i]
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
          n_loaded,
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

  # Aggregate short_visual counts per tumor
  if (length(visual_list) > 0) {
    visual_table <- bind_rows(visual_list) |>
      group_by(tumor_id, pattern, short_visual, R, L, prefix) |>
      summarise(count = sum(n), .groups = "drop")

    # Add cancer_type from sample_info
    if (!is.null(sample_info)) {
      visual_table$cancer_type <- sample_info$Cancer_Type[
        match(visual_table$tumor_id, sample_info$Patient)
      ]
    } else {
      visual_table$cancer_type <- NA_character_
    }

    # Compute R_intuitive: count how many times the deleted sequence
    # (between < and >) appears in the bracket context (between [ and ]).
    # Skip microhomology cases (curly braces inside <...>): set to NA.
    is_mh <- grepl("<\\{", visual_table$short_visual)
    del_seq <- sub(".*<([^>]+)>.*", "\\1", visual_table$short_visual)
    bracket_seq <- sub(".*\\[([^]]+)\\].*", "\\1", visual_table$short_visual)
    visual_table$R_intuitive <- mapply(function(d, b, mh) {
      if (mh) return(NA_integer_)
      if (nchar(d) == 0 || nchar(b) == 0) return(NA_integer_)
      m <- gregexpr(d, b, fixed = TRUE)[[1]]
      if (m[1] == -1L) 0L else length(m)
    }, del_seq, bracket_seq, is_mh, USE.NAMES = FALSE)

    # Microhomology length: number of characters between { and }
    visual_table$mh_length <- ifelse(
      grepl("\\{", visual_table$short_visual),
      nchar(sub(".*\\{([^}]+)\\}.*", "\\1", visual_table$short_visual)),
      NA_integer_
    )

    visual_table <- visual_table |>
      dplyr::select(
        tumor_id, cancer_type, pattern, short_visual,
        R, R_intuitive, L, mh_length, prefix, count
      ) |>
      arrange(tumor_id, pattern, desc(count))
  } else {
    visual_table <- tibble(
      tumor_id = character(),
      cancer_type = character(),
      pattern = character(),
      short_visual = character(),
      R = integer(),
      R_intuitive = integer(),
      L = integer(),
      mh_length = integer(),
      prefix = character(),
      count = integer()
    )
  }

  visual_table
}
