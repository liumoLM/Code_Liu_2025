# Helper functions for vignette.Rmd
# These functions separate computation from plotting for easier debugging

#" Format signature name with Greek letters
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

#' Create a fenced div with a style class
#'
#' @param txt Character: the text content for the div
#' @param style Character: the CSS class for the div (default ".callout-note")
#' @return Character string with fenced div markdown
fenced_div <- function(txt, style = ".callout-note") {
  paste0('\n\n::: {', style, '}\n', txt, '\n:::\n\n')
}

#' Output a fenced div to the document
#'
#' Convenience wrapper that combines cat() and fenced_div()
#'
#' @param txt Character: the text content for the div
#' @param style Character: the CSS class for the div (default ".callout-note")
#' @return NULL (called for side effect of outputting to document)
catfdiv <- function(txt, style = ".callout-note") {
  cat(fenced_div(txt, style))
}

#' Find signature-specific text file
#'
#' @param sig_id Character: the signature identifier
#' @return Character string with file contents, or NULL if file doesn't exist
find_sig_txt <- function(sig_id) {
  file_path <- file.path(data_dir, "per_sig_txt", glue::glue("{sig_id}.md"))
  if (!file.exists(file_path)) {
    return(NULL)
  }
  readLines(file_path, warn = FALSE) |> paste(collapse = "\n")
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

#' Generate markdown header text for a signature section
#'
#' @param sig_data List returned from compute_signature_data
#' @param use_html Logical: use HTML styling (for Quarto) or plain markdown
#' @return Character string with markdown/HTML text
generate_section_header <- function(sig_data, use_html = TRUE) {
  # Format signature name with Greek letters
  display_name <- format_signature_name(sig_data$ID89signature)
  firsttext = "corresponds to 83-type signature "

  if (use_html) {
    # Markdown heading for TOC entry, followed by styled content
    # The heading gets picked up by TOC, then we add styled div

    if (FALSE) {
      header <- paste0(
        '\n\n# ',
        display_name,
        ' {.signature-section}\n\n',
        '<div class="signature-header">',
        '<span class="signature-name">',
        display_name,
        '</span>',
        '</div>\n\n'
      )
    }

    header = ""

    if (sig_data$is_insdel15_16) {
      header <- paste0(
        header,
        '<div class="signature-subtitle">',
        '<strong>83-type signature:</strong> ',
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
          firsttext,
          sig_data$catalog,
          '<br><em>(Deletions and insertions in long poly T tracts were removed from the Indel83 spectrum)</em>',
          '</div>\n\n'
        )
      } else {
        header <- paste0(
          header,
          '<div class="signature-subtitle">',
          firsttext,
          sig_data$ID83signature,
          '</div>\n\n'
        )
      }

      if (sig_data$ID89signature == "InsDel_N") {
        header <- paste0(
          header,
          fenced_div(
            'InsDel_N is identical to InsDel_J in 89-Type representation, but different in 476-Type representation.',
            style = ".callout-note"
          )
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
  cosine476_text <- if (is.na(sig_data$cosine476)) {
    "N/A"
  } else {
    as.character(sig_data$cosine476)
  }

  df <- data.frame(
    `83-type` = sig_data$cosine83,
    `476-type` = cosine476_text,
    `89-type` = sig_data$cosine89,
    check.names = FALSE
  )

  if (use_html) {
    table_html <- knitr::kable(df, format = "html")
    # Add styling with internal gridlines and larger text
    table_output <- gsub(
      "<table>",
      '<table style="font-size: 0.75em; border-collapse: collapse; border: 1px solid black;"><style>table th, table td { border: 1px solid black; padding: 4px 8px; }</style>',
      table_html
    )
  } else {
    table_output <- knitr::kable(df, format = "markdown")
  }

  paste0(
    "\n\n",
    paste(table_output, collapse = "\n"),
    "\n\n---\n"
  )
}


#' Generate all plots for a signature and save to files
#'
#' @param sig_data List returned from compute_signature_data
#' @param ID89_signatures Data frame of ID89 signatures
#' @param ID89_catalogs Data frame of ID89 catalogs
#' @param ID83_signatures ICAMS catalog of ID83 signatures
#' @param ID83_catalogs ICAMS catalog of ID83 catalogs
#' @param ID83_catalogs_no_polyT ICAMS catalog with polyT removed
#' @param ID476_signatures Data frame of ID476 signatures
#' @param ID476_catalogs Data frame of ID476 catalogs
#' @param plot_dir Directory to save plots
#' @param plot476_base_size Base font size for 476 plots
#' @param plot476_label_size Label size for 476 plots
#' @param plot476_simplify_labels Whether to simplify labels
#' @return List with paths to all generated plot files
generate_plots_to_files <- function(
  sig_data,
  ID89_signatures,
  ID89_catalogs,
  ID83_signatures,
  ID83_catalogs,
  ID83_catalogs_no_polyT,
  ID476_signatures,
  ID476_catalogs,
  plot_dir,
  plot476_base_size = 20,
  plot476_label_size = 3,
  plot476_simplify_labels = FALSE
) {
  # Create safe filename prefix from signature name
  safe_name <- gsub("[^a-zA-Z0-9_]", "_", sig_data$ID89signature)

  paths <- list(
    id89_sig = file.path(plot_dir, paste0(safe_name, "_id89_sig.png")),
    id89_catalog = file.path(plot_dir, paste0(safe_name, "_id89_catalog.png")),
    id89_reconstructed = file.path(
      plot_dir,
      paste0(safe_name, "_id89_reconstructed.png")
    ),
    id89_diff = file.path(plot_dir, paste0(safe_name, "_id89_diff.png")),
    id476_sig = file.path(plot_dir, paste0(safe_name, "_id476_sig.png")),
    id476_catalog = file.path(
      plot_dir,
      paste0(safe_name, "_id476_catalog.png")
    ),
    id83_sig = file.path(plot_dir, paste0(safe_name, "_id83_sig.png")),
    id83_catalog = file.path(plot_dir, paste0(safe_name, "_id83_catalog.png"))
  )

  # Helper to save ggplot
  save_ggplot <- function(p, path, width = 19, height = 3) {
    ggplot2::ggsave(
      path,
      p,
      width = width,
      height = height,
      dpi = 150,
      bg = "white"
    )
  }

  # Helper to save base R plot (for ICAMS)
  save_base_plot <- function(expr, path, width = 10, height = 4) {
    png(
      path,
      width = width,
      height = height,
      units = "in",
      res = 150,
      bg = "white"
    )
    tryCatch(
      {
        force(expr)
      },
      finally = {
        dev.off()
      }
    )
  }

  p89 <- function(catalog, plot_title, setyaxis = NULL) {
    mSigPlot::plot_89(
      catalog,
      plot_title = plot_title,
      text_size = 5,
      base_size = getp('basesize89'),
      setyaxis = setyaxis
    )
  }

  save89 = function(myplot, path) {
    save_ggplot(myplot, path, width = getp('w89'), height = getp('h89'))
  }

  # ID89 Plot 1: Signature
  p1 <- p89(
    ID89_signatures[, sig_data$ID89signature, drop = FALSE],
    plot_title = sig_data$ID89signature
  )
  save89(p1, paths$id89_sig)

  # ID89 Plot 2: Catalog
  catalogtoplot = ID89_catalogs[, sig_data$catalog, drop = FALSE]
  ymax = max(catalogtoplot)
  p2 <- p89(
    catalogtoplot,
    plot_title = paste0(
      "Spectrum A: Mutational spectrum of ",
      sig_data$catalog,
      " | Cosine Similarity to ",
      sig_data$ID89signature,
      " = ",
      sig_data$cosine89
    ),
    setyaxis = ymax
  )
  save89(p2, paths$id89_catalog)

  # ID89 Plots 3 & 4: Decomposition (only for non-InsDel15/16)
  if (!sig_data$is_insdel15_16 && !is.null(sig_data$reconstructed_catalog)) {
    p3 <- p89(
      sig_data$reconstructed_catalog,
      plot_title = paste0(
        "Spectrum B: Partial mutational spectrum of ",
        sig_data$catalog,
        " due to ",
        sig_data$ID89signature
      ),
      setyaxis = ymax
    )
    save89(p3, paths$id89_reconstructed)

    p4 <- p89(
      sig_data$diff_catalog,
      plot_title = paste0(
        "Remaining mutations in ",
        sig_data$catalog,
        " not due to ",
        sig_data$ID89signature,
        " (A minus B) | Cosine similarity to ",
        sig_data$ID89signature,
        " = ",
        sig_data$cosine89_diff
      ),
      setyaxis = ymax
    )
    save89(p4, paths$id89_diff)
  } else {
    paths$id89_reconstructed <- NULL
    paths$id89_diff <- NULL
  }

  # ID476 plots
  p476 <- function(catalog, plot_title) {
    mSigPlot::plot_476(
      catalog,
      plot_title = plot_title,
      text_size = 5,
      label_size = plot476_label_size,
      num_labels = 5,
      base_size = plot476_base_size,
      simplify_labels = plot476_simplify_labels
    )
  }

  save476 = function(myplot, path) {
    save_ggplot(myplot, path, width = getp('w476'), height = getp('h476'))
  }

  if (sig_data$has_476_signature) {
    p5 <- p476(
      ID476_signatures[, sig_data$ID89signature],
      plot_title = paste0("476-Type Representation of ", sig_data$ID89signature)
    )
    save476(p5, paths$id476_sig)

    p6 <- p476(
      ID476_catalogs[, sig_data$catalog],
      plot_title = paste0("476-Type Spectrum of ", sig_data$catalog)
    )
    save476(p6, paths$id476_catalog)
  } else {
    p5 <- p476(
      ID476_catalogs[, sig_data$catalog],
      plot_title = paste0(
        "476-Type Representation of the Supporting Genome ",
        sig_data$catalog
      )
    )
    save476(p5, paths$id476_sig)
    paths$id476_catalog <- NULL
  }

  # ID83 signature (only if exists)

  p83 <- function(catalog, plot_title = NULL) {
    mSigPlot::plot_83(
      catalog,
      plot_title = plot_title,
      text_size = 5,
      base_size = getp('basesize83')
    )
  }
  save83 = function(myplot, path) {
    save_ggplot(myplot, path, width = getp('w83'), height = getp('h83'))
  }

  if (sig_data$has_83_signature) {
    ptmp = p83(ID83_signatures[,
      sig_data$ID83signature,
      drop = FALSE
    ])
    save83(ptmp, paths$id83_sig)
  } else {
    paths$id83_sig <- NULL
  }

  # ID83 spectrum catalog, always shown
  if (sig_data$is_polyT_removed) {
    cat83touse = ID83_catalogs_no_polyT
  } else {
    cat83touse = ID83_catalogs
  }
  ptmp <- p83(cat83touse[, sig_data$catalog, drop = FALSE])
  save83(ptmp, paths$id83_catalog)

  return(paths)
}


#' Generate all plots in parallel
#'
#' @param all_signature_data List of signature data from compute_signature_data
#' @param ... Additional arguments passed to generate_plots_to_files
#' @param n_workers Number of parallel workers (default 10)
#' @return List of plot paths for each signature
generate_all_plots_parallel <- function(
  all_signature_data,
  ID89_signatures,
  ID89_catalogs,
  ID83_signatures,
  ID83_catalogs,
  ID83_catalogs_no_polyT,
  ID476_signatures,
  ID476_catalogs,
  plot_dir,
  plot476_base_size = 20,
  plot476_label_size = 3,
  plot476_simplify_labels = FALSE,
  n_workers = 14
) {
  # Create plot directory
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

  # Set up parallel backend
  future::plan(future::multisession, workers = n_workers)

  # Generate plots in parallel
  all_paths <- furrr::future_map(
    all_signature_data,
    function(sig_data) {
      generate_plots_to_files(
        sig_data = sig_data,
        ID89_signatures = ID89_signatures,
        ID89_catalogs = ID89_catalogs,
        ID83_signatures = ID83_signatures,
        ID83_catalogs = ID83_catalogs,
        ID83_catalogs_no_polyT = ID83_catalogs_no_polyT,
        ID476_signatures = ID476_signatures,
        ID476_catalogs = ID476_catalogs,
        plot_dir = plot_dir,
        plot476_base_size = plot476_base_size,
        plot476_label_size = plot476_label_size,
        plot476_simplify_labels = plot476_simplify_labels
      )
    },
    .options = furrr::furrr_options(
      seed = TRUE,
      packages = c("ggplot2", "ICAMS", "mSigPlot", "indelsig.tools.lib")
    ),
    .progress = TRUE
  )

  # Reset to sequential
  future::plan(future::sequential)

  names(all_paths) <- names(all_signature_data)
  return(all_paths)
}


#' Check if plot cache is valid
#'
#' Compares hash of source data files against stored hash.
#' Returns TRUE if cache is valid (no regeneration needed).
#'
#' @param data_dir Directory containing source data files
#' @param plot_dir Directory where plots are stored
#' @param cache_file Name of the cache hash file
#' @return Logical: TRUE if cache is valid, FALSE if regeneration needed
check_plot_cache <- function(
  data_dir,
  plot_dir,
  cache_file = "plot_cache_hash.rds"
) {
  cache_path <- file.path(plot_dir, cache_file)

  # Source files that plots depend on
  source_files <- c(
    "Liu_et_al_final_89_type_signatures.tsv",
    "Liu_et_al_89_type_spectra.tsv",
    "Liu_et_al_final_83_type_signatures.tsv",
    "Liu_et_al_83_type_spectra.tsv",
    "Liu_et_al_final_476_type_signatures.tsv",
    "Liu_et_al_476_type_spectra.tsv",
    "89type_to_83type_connection.tsv"
  )

  # Check if plot directory exists
  if (!dir.exists(plot_dir)) {
    return(FALSE)
  }

  # Compute current hash from file modification times
  file_paths <- file.path(data_dir, source_files)
  if (!all(file.exists(file_paths))) {
    return(FALSE)
  }

  current_hash <- digest::digest(
    sapply(file_paths, file.mtime)
  )

  # Check if cache exists and matches
  if (file.exists(cache_path)) {
    stored_hash <- readRDS(cache_path)
    if (identical(stored_hash, current_hash)) {
      return(TRUE) # Cache is valid
    }
  }

  return(FALSE) # Cache invalid or missing
}


#' Save cache hash after generating plots
#'
#' @param data_dir Directory containing source data files
#' @param plot_dir Directory where plots are stored
#' @param cache_file Name of the cache hash file
save_plot_cache <- function(
  data_dir,
  plot_dir,
  cache_file = "plot_cache_hash.rds"
) {
  source_files <- c(
    "Liu_et_al_final_89_type_signatures.tsv",
    "Liu_et_al_89_type_spectra.tsv",
    "Liu_et_al_final_83_type_signatures.tsv",
    "Liu_et_al_83_type_spectra.tsv",
    "Liu_et_al_final_476_type_signatures.tsv",
    "Liu_et_al_476_type_spectra.tsv",
    "89type_to_83type_connection.tsv"
  )

  current_hash <- digest::digest(
    sapply(file.path(data_dir, source_files), file.mtime)
  )

  saveRDS(current_hash, file.path(plot_dir, cache_file))
}


#' Reconstruct plot paths from existing files
#'
#' Used when cache is valid to get paths without regenerating plots.
#'
#' @param signature_names Vector of signature names
#' @param plot_dir Directory where plots are stored
#' @return List of plot paths organized by signature name
reconstruct_plot_paths <- function(signature_names, plot_dir) {
  all_paths <- lapply(signature_names, function(sig_name) {
    safe_name <- gsub("[^a-zA-Z0-9_]", "_", sig_name)

    paths <- list(
      id89_sig = file.path(plot_dir, paste0(safe_name, "_id89_sig.png")),
      id89_catalog = file.path(
        plot_dir,
        paste0(safe_name, "_id89_catalog.png")
      ),
      id89_reconstructed = file.path(
        plot_dir,
        paste0(safe_name, "_id89_reconstructed.png")
      ),
      id89_diff = file.path(plot_dir, paste0(safe_name, "_id89_diff.png")),
      id476_sig = file.path(plot_dir, paste0(safe_name, "_id476_sig.png")),
      id476_catalog = file.path(
        plot_dir,
        paste0(safe_name, "_id476_catalog.png")
      ),
      id83_sig = file.path(plot_dir, paste0(safe_name, "_id83_sig.png")),
      id83_catalog = file.path(plot_dir, paste0(safe_name, "_id83_catalog.png"))
    )

    # Set to NULL if file doesn't exist
    paths <- lapply(paths, function(p) {
      if (file.exists(p)) p else NULL
    })

    return(paths)
  })

  names(all_paths) <- signature_names
  return(all_paths)
}
