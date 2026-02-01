#' Create overview kable table with hyperlinks to signature sections
#'
#' Generates a kable table for the Overview section with:
#' - Signature IDs hyperlinked to their corresponding sections
#' - User-friendly column names
#' - Centered text, smaller font
#' - Scrollable with fixed header
#'
#' @param sig_table Data frame: the signature summary table (sig_table2)
#' @return kableExtra styled table
create_overview_table <- function(sig_table) {
  source("table_1_col_name_mapping.R")
  source("vhelpers.R")

  df <- sig_table

  # Remove is_polyT_removed column (not needed in display)
  df$is_polyT_removed <- NULL

  # Remove the last column
  df <- df[, -ncol(df)]

  # Clean up exemplar_id: remove everything up to and including '::'
  df$exemplar_id <- sub(".*::", "", df$exemplar_id)

  # Clean up best_match_jin: remove 'jin' prefix
  df$best_match_jin <- sub("^jin", "", df$best_match_jin)

  # Round cos_v_koh first (before any HTML wrapping)
  df$cos_v_koh <- round(df$cos_v_koh, 4)

  # Handle best_match_koh duplicates: mark best match with asterisk, gray out non-best
  koh_values <- df$best_match_koh
  non_na_idx <- which(!is.na(koh_values))
  rows_to_gray <- c()

  unique_koh <- unique(koh_values[non_na_idx])
  for (koh_val in unique_koh) {
    matching_rows <- which(koh_values == koh_val)
    if (length(matching_rows) > 1) {
      # Duplicate: find row with highest cos_v_koh
      cos_vals <- df$cos_v_koh[matching_rows]
      best_idx <- matching_rows[which.max(cos_vals)]
      non_best_idx <- setdiff(matching_rows, best_idx)
      # Add asterisk to best match
      df$best_match_koh[best_idx] <- paste0(df$best_match_koh[best_idx], "*")
      # Track non-best rows for graying
      rows_to_gray <- c(rows_to_gray, non_best_idx)
    }
  }

  # Gray out non-best duplicate rows for best_match_koh and cos_v_koh
  if (length(rows_to_gray) > 0) {
    df$best_match_koh[rows_to_gray] <- paste0(
      '<span style="color: #B0B0B0;">', df$best_match_koh[rows_to_gray], '</span>'
    )
    df$cos_v_koh[rows_to_gray] <- paste0(
      '<span style="color: #B0B0B0;">', df$cos_v_koh[rows_to_gray], '</span>'
    )
  }

  # Add hyperlinks to signature IDs
  df$signature_id <- sapply(df$signature_id, make_sig_hyperlink)

  # Rename columns using the mapping
  names(df) <- table_1_col_name_mapping(names(df))

  # Round numeric columns to 4 decimal places
  numeric_cols <- sapply(df, is.numeric)
  df[numeric_cols] <- lapply(df[numeric_cols], function(x) round(x, 4))

  # Create kable with HTML formatting
  tbl <- knitr::kable(df, format = "html", escape = FALSE, align = 'c') |>
    kableExtra::kable_styling(
      bootstrap_options = c("striped", "hover", "condensed"),
      full_width = TRUE,
      fixed_thead = TRUE
    ) |>
    kableExtra::scroll_box(width = "100%", height = "600px")

  # Wrap in div with smaller font size
  htmltools::div(style = "font-size: 70%;", htmltools::HTML(tbl))
}
