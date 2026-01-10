# Helper functions for vignette.Rmd
# These functions separate computation from plotting for easier debugging

#' Format signature name with Greek letters
#'
#' Replaces _alpha with α and _beta with β in signature names
#'
#' @param name Character: the signature name
#' @return Character: formatted name with Greek letters
format_signature_name <- function(name) {
  name <- gsub("_alpha", "α", name)
  name <- gsub("_beta", "β", name)
  name
}

#' Compute cosine similarities for a signature-catalog pair
#'
#' @param ID89signature Character: the ID89 signature name
#' @param catalog Character: the catalog/sample name
#' @param ID89_signatures Data frame of ID89 signatures
#' @param ID89_catalogs Data frame of ID89 catalogs
#' @param ID83_signatures ICAMS catalog of ID83 signatures
#' @param ID83_catalogs ICAMS catalog of ID83 catalogs
#' @param ID83_catalogs_no_polyT ICAMS catalog with polyT removed
#' @param ID476_signatures Data frame of ID476 signatures
#' @param ID476_catalogs Data frame of ID476 catalogs
#' @param ID83signature Character: the corresponding ID83 signature name
#' @param assignment_matrix Data frame of signature assignments
#' @return List with cosine similarities and intermediate data
compute_signature_data <- function(
  ID89signature,
  catalog,
  ID83signature,
  ID89_signatures,
  ID89_catalogs,
  ID83_signatures,
  ID83_catalogs,
  ID83_catalogs_no_polyT,
  ID476_signatures,
  ID476_catalogs,
  assignment_matrix
) {
  message("ID89signature = ", ID89signature)
  result <- list(
    ID89signature = ID89signature,
    catalog = catalog,
    ID83signature = ID83signature,
    is_insdel15_16 = ID89signature %in% c("InsDel15", "InsDel16"),
    is_polyT_removed = ID83signature %in%
      c("C_ID7", "ID_J", "C_ID10", "ID_N", "ID_O"),
    has_476_signature = ID89signature %in% colnames(ID476_signatures),
    has_83_signature = ID83signature %in% colnames(ID83_signatures)
  )

  # Compute cosine89 (raw catalog vs signature)
  result$cosine89 <- round(
    lsa::cosine(
      as.numeric(ID89_signatures[, ID89signature]),
      as.numeric(ID89_catalogs[, catalog])
    ),
    3
  )

  # For non-InsDel15/16, compute the decomposition

  if (!result$is_insdel15_16) {
    # Get assignment for this catalog
    assignment <- assignment_matrix[, catalog, drop = FALSE]

    # Zero out the current signature to get "other signatures" contribution
    assignment_others <- assignment
    assignment_others[
      which(row.names(assignment_matrix) == ID89signature),
    ] <- 0

    # Reconstruct catalog without this signature
    # Using columns 1:44 and 46 (excluding column 45)
    result$reconstructed_catalog <- as.matrix(ID89_signatures[, c(
      1:44,
      46
    )]) %*%
      as.matrix(assignment_others)

    # Difference = mutations attributed to this signature
    result$diff_catalog <- ID89_catalogs[, catalog, drop = FALSE] -
      result$reconstructed_catalog
    result$diff_catalog[result$diff_catalog < 0] <- 0

    # Cosine of diff vs signature
    result$cosine89_diff <- round(
      lsa::cosine(
        as.numeric(ID89_signatures[, ID89signature]),
        as.numeric(as.matrix(result$diff_catalog))
      ),
      3
    )
  } else {
    result$cosine89_diff <- NA
    result$reconstructed_catalog <- NULL
    result$diff_catalog <- NULL
  }

  # Compute cosine476
  if (result$has_476_signature) {
    result$cosine476 <- round(
      lsa::cosine(
        as.numeric(ID476_signatures[, ID89signature]),
        as.numeric(ID476_catalogs[, catalog])
      ),
      3
    )
  } else {
    result$cosine476 <- NA
  }

  # Compute cosine83
  if (result$is_polyT_removed) {
    result$cosine83 <- round(
      lsa::cosine(
        as.numeric(ID83_signatures[, ID83signature]),
        as.numeric(ID83_catalogs_no_polyT[, catalog])
      ),
      3
    )
  } else {
    if (ID83signature %in% colnames(ID83_signatures)) {
      result$cosine83 <- round(
        lsa::cosine(
          as.numeric(ID83_signatures[, ID83signature]),
          as.numeric(ID83_catalogs[, catalog])
        ),
        3
      )
    } else {
      result$cosine83 <- 0
    }
  }

  return(result)
}


#' Create plots for ID89 signature visualization
#'
#' @param sig_data List returned from compute_signature_data
#' @param ID89_signatures Data frame of ID89 signatures
#' @param ID89_catalogs Data frame of ID89 catalogs
#' @param plot89_height Unit for plot height
#' @return List of ggplot objects
create_id89_plots <- function(
  sig_data,
  ID89_signatures,
  ID89_catalogs,
  plot89_height
) {
  plots <- list()

  # Plot 1: The signature itself
  plots$p1 <- plot_89(
    ID89_signatures[, sig_data$ID89signature, drop = FALSE],
    text_size = 5,
    plot_title = sig_data$ID89signature,
    ylabel = "Props"
  )

  # Plot 2: The catalog spectrum
  plots$p2 <- plot_89(
    ID89_catalogs[, sig_data$catalog, drop = FALSE],
    text_size = 5,
    plot_title = paste0(
      "Spectrum A: Mutational spectrum of ",
      sig_data$catalog,
      " | Cosine Similarity to ",
      sig_data$ID89signature,
      " = ",
      sig_data$cosine89
    )
  )

  # For non-InsDel15/16, add decomposition plots
  if (!sig_data$is_insdel15_16 && !is.null(sig_data$reconstructed_catalog)) {
    plots$p3 <- plot_89(
      sig_data$reconstructed_catalog,
      text_size = 5,
      plot_title = paste0(
        "Spectrum B: Partial mutational spectrum of ",
        sig_data$catalog,
        " not using ",
        sig_data$ID89signature
      ),
      setyaxis = max(sig_data$diff_catalog)
    )

    plots$p4 <- plot_89(
      sig_data$diff_catalog,
      text_size = 5,
      plot_title = paste0(
        "Mutations due to ",
        sig_data$ID89signature,
        " (A minus B) | Cosine similarity to ",
        sig_data$ID89signature,
        " = ",
        sig_data$cosine89_diff
      )
    )
  }

  return(plots)
}


#' Create plots for ID476 signature visualization
#'
#' @param sig_data List returned from compute_signature_data
#' @param ID476_signatures Data frame of ID476 signatures
#' @param ID476_catalogs Data frame of ID476 catalogs
#' @param plot476_base_size Base font size for plot
#' @param plot476_simplify_labels Logical for label simplification
#' @return List of ggplot objects
create_id476_plots <- function(
  sig_data,
  ID476_signatures,
  ID476_catalogs,
  plot476_base_size,
  plot476_simplify_labels
) {
  p476 = function(catalog, plot_title) {
    plot_476(
      catalog,
      plot_title = plot_title,
      text_size = 5,
      label_size = plot476_label_size,
      num_labels = 5,
      base_size = plot476_base_size,
      simplify_labels = plot476_simplify_labels
    )
  }

  plots <- list()

  if (sig_data$has_476_signature) {
    plots$p5 <- p476(
      ID476_signatures[, sig_data$ID89signature],
      plot_title = paste0(
        "476-Type Representation of ",
        sig_data$ID89signature
      )
    )

    plots$p6 <- p476(
      ID476_catalogs[, sig_data$catalog],
      plot_title = paste0(
        "476-Type Spectrum of ",
        sig_data$catalog
      )
    )
  } else {
    # No extracted 476 signature, show catalog only
    plots$p5 <- p476(
      ID476_catalogs[, sig_data$catalog],
      plot_title = paste0(
        "476-Type Representation of the Supporting Genome ",
        sig_data$catalog
      )
    )
    plots$p6 <- NULL
  }

  return(plots)
}


#' Generate markdown header text for a signature section
#'
#' @param sig_data List returned from compute_signature_data
#' @param use_html Logical: use HTML styling (for Quarto) or plain markdown
#' @return Character string with markdown/HTML text
generate_section_header <- function(sig_data, use_html = TRUE) {
  # Format signature name with Greek letters
  display_name <- format_signature_name(sig_data$ID89signature)

  if (use_html) {
    # Markdown heading for TOC entry, followed by styled content
    # The heading gets picked up by TOC, then we add styled div
    header <- paste0(
      '\n\n### ',
      display_name,
      ' {.signature-section}\n\n',
      '<div class="signature-header">',
      '<span class="signature-name">',
      display_name,
      '</span>',
      '</div>\n\n'
    )

    if (sig_data$is_insdel15_16) {
      header <- paste0(
        header,
        '<div class="signature-subtitle">',
        '<strong>83-Type Signature:</strong> ',
        sig_data$ID83signature,
        '<br><strong>Note:</strong> The signature contributes all mutations of the example spectrum: ',
        sig_data$catalog,
        '</div>\n\n'
      )
    } else {
      if (sig_data$is_polyT_removed) {
        header <- paste0(
          header,
          '<div class="signature-subtitle">',
          '<strong>Supporting spectrum:</strong> ',
          sig_data$catalog,
          '<br><em>(Deletions and insertions in long poly T tracts were removed from the Indel83 spectrum)</em>',
          '</div>\n\n'
        )
      } else {
        header <- paste0(
          header,
          '<div class="signature-subtitle">',
          '<strong>83-type signature:</strong> ',
          sig_data$ID83signature,
          '<br><strong>Supporting spectrum:</strong> ',
          sig_data$catalog,
          '</div>\n\n'
        )
      }

      if (sig_data$ID89signature == "InsDel_N") {
        header <- paste0(
          header,
          '<div class="signature-note">',
          'InsDel_N is identical to InsDel_J in 89-Type representation, ',
          'but different in 476-Type representation.',
          '</div>\n\n'
        )
      }
    }
  } else {
    # Plain markdown fallback
    header <- paste0("\n\n### ", display_name, "\n")

    if (sig_data$is_insdel15_16) {
      header <- paste0(
        header,
        "\n83-type signature: ",
        sig_data$ID83signature,
        "\n\nThe signature contributes all mutations of the example spectrum: ",
        sig_data$catalog,
        "\n\n"
      )
    } else {
      if (sig_data$is_polyT_removed) {
        header <- paste0(
          header,
          "\n\n#### Supporting spectrum (deletions and insertions in long poly T tracts were removed from the Indel83 spectrum): ",
          sig_data$catalog,
          "\n\n"
        )
      } else {
        header <- paste0(
          header,
          "\n83-type signature: ",
          sig_data$ID83signature,
          "\n\nSupporting spectrum: ",
          sig_data$catalog,
          "\n\n"
        )
      }

      if (sig_data$ID89signature == "InsDel_N") {
        header <- paste0(
          header,
          "\nInsDel_N is identical to InsDel_J in 89-type representation, ",
          "but different in 476-type representation\n"
        )
      }
    }
  }

  return(header)
}


#' Generate markdown footer text with cosine summary
#'
#' @param sig_data List returned from compute_signature_data
#' @param use_html Logical: use HTML styling (for Quarto) or plain markdown
#' @return Character string with markdown/HTML text
generate_section_footer <- function(sig_data, use_html = TRUE) {
  cosine476_text <- if (is.na(sig_data$cosine476)) "N/A" else sig_data$cosine476

  if (use_html) {
    paste0(
      '\n\n<div class="cosine-summary">',
      '<h4>Cosine Similarities Summary</h4>',
      '<ul>',
      '<li><strong>83-type:</strong> <span class="cosine-value">',
      sig_data$cosine83,
      '</span></li>',
      '<li><strong>89-type:</strong> <span class="cosine-value">',
      sig_data$cosine89,
      '</span></li>',
      '<li><strong>476-type:</strong> <span class="cosine-value">',
      cosine476_text,
      '</span></li>',
      '</ul>',
      '</div>\n\n',
      '<hr class="section-divider">\n\n'
    )
  } else {
    paste0(
      "\n\n#### Cosine Similarities Summary\n",
      "- 83-type: ",
      sig_data$cosine83,
      "\n- 89-type: ",
      sig_data$cosine89,
      "\n- 476-type: ",
      cosine476_text,
      "\n\n---\n"
    )
  }
}
