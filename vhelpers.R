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
#' @param exemplar_id Character: the identifier of the supporting tumor
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
    exemplar_id,
    ID83signature,
    ID89_signatures,
    ID89_catalogs,
    ID83_signatures,
    ID83_catalogs,
    ID83_catalogs_no_polyT,
    ID476_signatures,
    ID476_catalogs,
    assignment_matrix,
    ID89_mapped_signatures = NULL,
    ID83_mapped_signatures = NULL
) {
  message("ID89signature = ", ID89signature)
  # Check if mapped 89-type signature exists (column name is {signature}_converted)
  mapped_col_name <- paste0(ID89signature, "_converted")
  has_mapped_sig <- !is.null(ID89_mapped_signatures) &&
    mapped_col_name %in% colnames(ID89_mapped_signatures)
  
  # Check if mapped 83-type signature exists
  # Only consider it if a 476-type signature was extracted
  has_476_sig <- ID89signature %in% colnames(ID476_signatures)
  has_83_mapped_sig <- has_476_sig &&
    !is.null(ID83_mapped_signatures) &&
    mapped_col_name %in% colnames(ID83_mapped_signatures)
  
  result <- list(
    ID89signature = ID89signature,
    catalog = exemplar_id,
    ID83signature = ID83signature,
    is_insdel15_16 = ID89signature %in% c("InsDel15", "InsDel16"),
    is_polyT_removed = ID83signature %in%
      c("C_ID7", "ID_J", "C_ID10", "ID_N", "ID_O"),
    has_476_signature = has_476_sig,
    has_83_signature = ID83signature %in% colnames(ID83_signatures),
    has_mapped_signature = has_mapped_sig,
    has_83_mapped_signature = has_83_mapped_sig
  )
  
  # Compute cosine89 (raw catalog vs signature)
  result$cosine89 <-
    lsa::cosine(
      as.numeric(ID89_signatures[, ID89signature]),
      as.numeric(ID89_catalogs[, exemplar_id])
    )
  
  # Compute cosine similarity between main signature and mapped signature
  if (has_mapped_sig) {
    result$cosine89_mapped <-
      lsa::cosine(
        as.numeric(ID89_signatures[, ID89signature]),
        as.numeric(ID89_mapped_signatures[, mapped_col_name])
      )
  } else {
    result$cosine89_mapped <- NA
  }
  
  # Assert that all signatures in assignment matrix exist in signature matrix
  missing_sigs <- setdiff(
    row.names(assignment_matrix),
    colnames(ID89_signatures)
  )
  if (length(missing_sigs) > 0) {
    stop(
      "dim assignment matrix ",
      dim(assignment_matrix),
      "Assignment matrix has signatures not in signature matrix: ",
      paste(missing_sigs, collapse = ", ")
    )
  }
  
  # For non-InsDel15/16, compute the decomposition
  if (!result$is_insdel15_16) {
    # Find signatures present in both the signature matrix and assignment matrix
    common_sigs <- intersect(
      colnames(ID89_signatures),
      row.names(assignment_matrix)
    )
    
    # Get assignment for this catalog, filtered to common signatures
    assignment <- assignment_matrix[common_sigs, exemplar_id, drop = FALSE]
    
    # There is no assignment for InsDel_N because it is identical to
    # InsDel_J
    sigid = ID89signature
    if (sigid == "InsDel_N") {
      sigid <- "InsDel_J"
    }
    # Zero out the current signature to get contribution by other signatures
    assignment_others <- assignment
    stopifnot(sigid %in% common_sigs)
    assignment_others[sigid, ] <- 0
    rm(sigid)
    
    # Reconstruct catalog without this signature
    # Use common signature names to ensure alignment
    result$residual_spectrum <- as.matrix(
      ID89_signatures[, common_sigs]
    ) %*%
      as.matrix(assignment_others)
    
    # Difference = mutations attributed to this signature
    result$target_sig_partial_spectrum <- ID89_catalogs[,
                                                        exemplar_id,
                                                        drop = FALSE
    ] -
      result$residual_spectrum
    result$target_sig_partial_spectrum[
      result$target_sig_partial_spectrum < 0
    ] <- 0
    
    # Cosine of diff vs signature
    result$cosine89_diff <-
      lsa::cosine(
        as.numeric(ID89_signatures[, ID89signature]),
        as.numeric(as.matrix(result$target_sig_partial_spectrum))
      )
  } else {
    result$cosine89_diff <- NA
    result$residual_spectrum <- NULL
    result$target_sig_partial_spectrum <- NULL
  }
  
  # Compute cosine476
  if (result$has_476_signature) {
    result$cosine476 <-
      lsa::cosine(
        as.numeric(ID476_signatures[, ID89signature]),
        as.numeric(ID476_catalogs[, exemplar_id])
      )
  } else {
    result$cosine476 <- NA
  }
  
  # Compute cosine83
  if (result$is_polyT_removed) {
    result$cosine83 <-
      lsa::cosine(
        as.numeric(ID83_signatures[, ID83signature]),
        as.numeric(ID83_catalogs_no_polyT[, exemplar_id])
      )
  } else {
    if (ID83signature %in% colnames(ID83_signatures)) {
      result$cosine83 <-
        lsa::cosine(
          as.numeric(ID83_signatures[, ID83signature]),
          as.numeric(ID83_catalogs[, exemplar_id])
        )
    } else {
      result$cosine83 <- 0
    }
  }
  
  # Compute cosine similarity between native 83-type and mapped 83-type signature
  if (has_83_mapped_sig && result$has_83_signature) {
    result$cosine83_mapped <-
      lsa::cosine(
        as.numeric(ID83_signatures[, ID83signature]),
        as.numeric(ID83_mapped_signatures[, mapped_col_name])
      )
  } else {
    result$cosine83_mapped <- NA
  }
  
  return(result)
}

#' Generate markdown footer text with cosine summary
#'
#' @param sig_data List returned from compute_signature_data
#' @return Character string with markdown text
generate_section_footer <- function(sig_data) {
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
  
  table_output <- knitr::kable(df)
  
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
    plot476_simplify_labels = FALSE,
    ID89_mapped_signatures = NULL,
    ID83_mapped_signatures = NULL
) {
  # Create safe filename prefix from signature name
  safe_name <- gsub("[^a-zA-Z0-9_]", "_", sig_data$ID89signature)
  
  paths <- list(
    id89_sig = file.path(plot_dir, paste0(safe_name, "_id89_sig.png")),
    id89_mapped = file.path(plot_dir, paste0(safe_name, "_id89_mapped.png")),
    id89_catalog = file.path(plot_dir, paste0(safe_name, "_id89_catalog.png")),
    id89_residual = file.path(
      plot_dir,
      paste0(safe_name, "_id89_residual.png")
    ),
    id89_target_sig_partial_spectrum = file.path(
      plot_dir,
      paste0(safe_name, "_id89_target_sig_partial_spectrum.png")
    ),
    id476_sig = file.path(plot_dir, paste0(safe_name, "_id476_sig.png")),
    id476_catalog = file.path(
      plot_dir,
      paste0(safe_name, "_id476_catalog.png")
    ),
    id83_sig = file.path(plot_dir, paste0(safe_name, "_id83_sig.png")),
    id83_mapped = file.path(plot_dir, paste0(safe_name, "_id83_mapped.png")),
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
  
  p89 <- function(catalog, plot_title, setyaxis = NULL) {
   p <- mSigPlot::plot_89(
      catalog,
      plot_title = plot_title,
      text_size = getp('textsize89'),
      base_size = getp('basesize89'),
      setyaxis = setyaxis
    )
    ##  🎯🎯🎯 Modification of Changing text 🎯🎯🎯  xueming
    for (i in seq_along(p$layers)) {
      layer <- p$layers[[i]]
      
      if (inherits(layer$geom, c("GeomText", "GeomLabel"))) {
        layer$aes_params$size   <- 6
        layer$aes_params$colour <- "black"  # or red for test
        p$layers[[i]] <- layer
      }
    }
    
    return(p)
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
  
  # ID89 Plot 1b: Mapped signature (from 476-type)
  if (sig_data$has_mapped_signature && !is.null(ID89_mapped_signatures)) {
    mapped_col_name <- paste0(sig_data$ID89signature, "_converted")
    p1b <- p89(
      ID89_mapped_signatures[, mapped_col_name, drop = FALSE],
      plot_title = paste0(
        sig_data$ID89signature,
        " converted from 476-type signature | cosine similarity to ",
        sig_data$ID89signature,
        " = ",
        format(sig_data$cosine89_mapped, digits = getp("cosine_digits"))
      )
    )
    save89(p1b, paths$id89_mapped)
  } else {
    paths$id89_mapped <- NULL
  }
  
  # ID89 Plot 2: Catalog
  catalogtoplot = ID89_catalogs[, sig_data$catalog, drop = FALSE]
  ymax = max(catalogtoplot)
  p2 <- p89(
    catalogtoplot,
    plot_title = paste0(
      "Spectrum A, from ",
      sig_data$catalog,
      " | cosine similarity to ",
      sig_data$ID89signature,
      " = ",
      format(sig_data$cosine89, digits = getp("cosine_digits"))
    ),
    setyaxis = ymax
  )
  save89(p2, paths$id89_catalog)
  
  # ID89 Plots 3 & 4: Decomposition (only for non-InsDel15/16)
  if (!sig_data$is_insdel15_16 && !is.null(sig_data$residual_spectrum)) {
    target_sig_title <- paste0(
      "Spectrum B: partial mutational spectrum of ",
      sig_data$catalog,
      " due to ",
      sig_data$ID89signature
    )
    
    residual_title <- paste0(
      "Remaining mutations in ",
      sig_data$catalog,
      " not due to ",
      sig_data$ID89signature,
      " (A minus B) " #| Cosine similarity to ",
      # sig_data$ID89signature,
      # " = ",
      # format(sig_data$cosine89_diff, digits = 4)
    )

    p3 <- p89(
      sig_data$residual_spectrum,
      plot_title = residual_title,
      setyaxis = ymax
    )
    save89(p3, paths$id89_residual)
    
    p4 <- p89(
      sig_data$target_sig_partial_spectrum,
      plot_title = target_sig_title,
      setyaxis = ymax
    )
    save89(p4, paths$id89_target_sig_partial_spectrum)
} else {
  paths$id89_residual <- NULL
  paths$id89_target_sig_partial_spectrum <- NULL
}
  
  # ID476 plots
  p476 <- function(catalog, plot_title) {
    p <- mSigPlot::plot_476(
      catalog,
      plot_title = plot_title,
      text_size = 5,
      label_size = plot476_label_size,
      num_labels = 5,
      base_size = plot476_base_size,
      simplify_labels = plot476_simplify_labels
    )
    
    
    # 🎯🎯🎯Repairing text in drawn ggplot objects  🎯🎯🎯xueming
    
    # 1. Repair master data
    if (!is.null(p$data)) {
      for (col in names(p$data)) {
        if (is.character(p$data[[col]])) {
          p$data[[col]] <- gsub("22bp", "2bp", p$data[[col]], fixed = TRUE)
        } else if (is.factor(p$data[[col]])) {
          levels(p$data[[col]]) <- gsub("22bp", "2bp", levels(p$data[[col]]), fixed = TRUE)
        }
      }
    }
    
    # 2. Repair layer data
    for (i in seq_along(p$layers)) {
      if (!is.null(p$layers[[i]]$data) && is.data.frame(p$layers[[i]]$data)) {
        for (col in names(p$layers[[i]]$data)) {
          if (is.character(p$layers[[i]]$data[[col]])) {
            p$layers[[i]]$data[[col]] <- gsub("22bp", "2bp", p$layers[[i]]$data[[col]], fixed = TRUE)
          } else if (is.factor(p$layers[[i]]$data[[col]])) {
            levels(p$layers[[i]]$data[[col]]) <- gsub("22bp", "2bp", levels(p$layers[[i]]$data[[col]]), fixed = TRUE)
          }
        }
      }
    }
    
    # 3. Fix the scale on the x and y axes.
    for (i in seq_along(p$scales$scales)) {
      if (!is.null(p$scales$scales[[i]]$breaks)) {
        if (is.character(p$scales$scales[[i]]$breaks)) {
          p$scales$scales[[i]]$breaks <- gsub("22bp", "2bp", p$scales$scales[[i]]$breaks, fixed = TRUE)
        }
      }
      if (!is.null(p$scales$scales[[i]]$labels)) {
        if (is.character(p$scales$scales[[i]]$labels) || is.function(p$scales$scales[[i]]$labels)) {
          if (is.character(p$scales$scales[[i]]$labels)) {
            p$scales$scales[[i]]$labels <- gsub("22bp", "2bp", p$scales$scales[[i]]$labels, fixed = TRUE)
          }
        }
      }
    }
    
    return(p)
  }
  
  save476 = function(myplot, path) {
    save_ggplot(myplot, path, width = getp('w476'), height = getp('h476'))
  }
  
  if (sig_data$has_476_signature) {
    p5 <- p476(
      ID476_signatures[, sig_data$ID89signature],
      plot_title = paste0(
        "Extracted 476-type signature corresponding to ",
        sig_data$ID89signature
      )
    )
    save476(p5, paths$id476_sig)
    
    p6 <- p476(
      ID476_catalogs[, sig_data$catalog],
      plot_title = ""
    )
    save476(p6, paths$id476_catalog)
  } else {
    p5 <- p476(
      ID476_catalogs[, sig_data$catalog],
      plot_title = paste0(
        "476-type spectrum of the supporting tumor ",
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
      text_size = getp('textsize83'),
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
  
  # ID83 mapped signature (from 476-type)
  if (sig_data$has_83_mapped_signature && !is.null(ID83_mapped_signatures)) {
    mapped_col_name <- paste0(sig_data$ID89signature, "_converted")
    cosine_text <- if (!is.na(sig_data$cosine83_mapped)) {
      paste0(
        " | cosine to native = ",
        format(sig_data$cosine83_mapped, digits = getp("cosine_digits"))
      )
    } else {
      ""
    }
    ptmp <- p83(
      ID83_mapped_signatures[, mapped_col_name, drop = FALSE],
      plot_title = ""
    )
    save83(ptmp, paths$id83_mapped)
  } else {
    paths$id83_mapped <- NULL
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
    ID89_mapped_signatures = NULL,
    ID83_mapped_signatures = NULL,
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
        plot476_simplify_labels = plot476_simplify_labels,
        ID89_mapped_signatures = ID89_mapped_signatures,
        ID83_mapped_signatures = ID83_mapped_signatures
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
  
  # Source files in data_dir that plots depend on
  data_files <- c(
    "Liu_et_al_final_89_type_signatures.tsv",
    "Liu_et_al_89_type_spectra.tsv",
    "Liu_et_al_final_83_type_signatures.tsv",
    "Liu_et_al_83_type_spectra.tsv",
    "Liu_et_al_final_476_type_signatures.tsv",
    "Liu_et_al_476_type_spectra.tsv",
    "89type_to_83type_connection.tsv"
  )
  
  # Files in vignette directory that affect plotting
  vignette_files <- c(
    "ppar.R",
    "89_mapped_from_476.tsv",
    "83_mapped_from_476.tsv"
  )
  
  # Check if plot directory exists
  if (!dir.exists(plot_dir)) {
    return(FALSE)
  }
  
  # Compute current hash from file modification times
  data_paths <- file.path(data_dir, data_files)
  if (!all(file.exists(data_paths))) {
    return(FALSE)
  }
  
  # Check vignette files (optional - skip missing files)
  vignette_paths <- vignette_files[file.exists(vignette_files)]
  
  all_paths <- c(data_paths, vignette_paths)
  current_hash <- digest::digest(
    sapply(all_paths, file.mtime)
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
  # Source files in data_dir
  data_files <- c(
    "Liu_et_al_final_89_type_signatures.tsv",
    "Liu_et_al_89_type_spectra.tsv",
    "Liu_et_al_final_83_type_signatures.tsv",
    "Liu_et_al_83_type_spectra.tsv",
    "Liu_et_al_final_476_type_signatures.tsv",
    "Liu_et_al_476_type_spectra.tsv",
    "89type_to_83type_connection.tsv"
  )
  
  # Files in vignette directory that affect plotting
  vignette_files <- c(
    "ppar.R",
    "89_mapped_from_476.tsv",
    "83_mapped_from_476.tsv"
  )
  
  data_paths <- file.path(data_dir, data_files)
  vignette_paths <- vignette_files[file.exists(vignette_files)]
  
  all_paths <- c(data_paths, vignette_paths)
  current_hash <- digest::digest(
    sapply(all_paths, file.mtime)
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
      id89_mapped = file.path(plot_dir, paste0(safe_name, "_id89_mapped.png")),
      id89_catalog = file.path(
        plot_dir,
        paste0(safe_name, "_id89_catalog.png")
      ),
      id89_residual = file.path(
        plot_dir,
        paste0(safe_name, "_id89_residual.png")
      ),
      id89_target_sig_partial_spectrum = file.path(
        plot_dir,
        paste0(safe_name, "_id89_target_sig_partial_spectrum.png")
      ),
      id476_sig = file.path(plot_dir, paste0(safe_name, "_id476_sig.png")),
      id476_catalog = file.path(
        plot_dir,
        paste0(safe_name, "_id476_catalog.png")
      ),
      id83_sig = file.path(plot_dir, paste0(safe_name, "_id83_sig.png")),
      id83_mapped = file.path(plot_dir, paste0(safe_name, "_id83_mapped.png")),
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
