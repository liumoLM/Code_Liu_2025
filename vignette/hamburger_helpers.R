# Helper functions for interactive hamburger plots
# These functions prepare data for plotly-based interactive visualization

#' Get hamburger plot data for a signature
#'
#' Extracts and prepares data for creating an interactive hamburger plot
#' showing mutation distribution across cancer types.
#'
#' @param signature_name Character: name of the signature (row name in matrix).
#' @param assignment_matrix Data frame with signatures as rows, samples as columns.
#'   Column names should be in format "CancerType::SampleID".
#' @param min_samples Integer: minimum samples per cancer type to include.
#' @param exclude_zero Logical: whether to exclude zero-mutation samples.
#' @return Data frame with columns: sample_id, cancer_type, mutations, x_pos,
#'   x_plot, median_mutations, and ordering information.
get_hamburger_data <- function(
  signature_name,
  assignment_matrix,
  min_samples = 3,
  exclude_zero = TRUE
) {
  if (!signature_name %in% rownames(assignment_matrix)) {
    stop("Signature '", signature_name, "' not found in assignment matrix")
  }

  sig_values <- as.numeric(assignment_matrix[signature_name, ])
  sample_names <- colnames(assignment_matrix)
  cancer_types <- sub("::.*", "", sample_names)

  # Create data frame with all samples

df <- data.frame(
    sample_id = sample_names,
    cancer_type = cancer_types,
    mutations = sig_values,
    stringsAsFactors = FALSE
  )

  # Calculate sample counts before filtering
  total_by_type <- as.data.frame(table(df$cancer_type))
  colnames(total_by_type) <- c("cancer_type", "total_samples")

  nonzero_by_type <- as.data.frame(table(df$cancer_type[df$mutations > 0]))
  colnames(nonzero_by_type) <- c("cancer_type", "nonzero_samples")

  sample_counts <- merge(total_by_type, nonzero_by_type, by = "cancer_type", all.x = TRUE)
  sample_counts$nonzero_samples[is.na(sample_counts$nonzero_samples)] <- 0

  # Filter zeros if requested
  if (exclude_zero) {
    df <- df[df$mutations > 0, ]
  }

  # Filter cancer types with insufficient samples
  type_counts <- table(df$cancer_type)
  valid_types <- names(type_counts)[type_counts >= min_samples]
  df <- df[df$cancer_type %in% valid_types, ]

  if (nrow(df) == 0) {
    return(NULL)
  }

  # Calculate medians for each cancer type
  medians <- aggregate(mutations ~ cancer_type, data = df, FUN = median)
  colnames(medians) <- c("cancer_type", "median_mutations")

  # Order by median
  medians <- medians[order(medians$median_mutations), ]
  type_order <- medians$cancer_type

  # Create ordered factor
  df$cancer_type <- factor(df$cancer_type, levels = type_order)
  df$x_pos <- as.numeric(df$cancer_type)

  # Sort points within each cancer type by mutation value
  df <- df[order(df$cancer_type, df$mutations), ]

  # Calculate within-group positions for spreading points
  df$within_rank <- ave(
    seq_len(nrow(df)),
    df$cancer_type,
    FUN = function(x) seq_along(x)
  )
  df$group_size <- ave(
    seq_len(nrow(df)),
    df$cancer_type,
    FUN = length
  )

  # Spread points within each cancer type column
  df$x_plot <- df$x_pos +
    0.8 * (df$within_rank - 1) / pmax(df$group_size - 1, 1) - 0.4
  df$x_plot[df$group_size == 1] <- df$x_pos[df$group_size == 1]

  # Add median info
  df <- merge(df, medians, by = "cancer_type")

  # Add sample counts
  sample_counts <- sample_counts[sample_counts$cancer_type %in% type_order, ]
  sample_counts$cancer_type <- factor(sample_counts$cancer_type, levels = type_order)
  df <- merge(df, sample_counts, by = "cancer_type")

  # Restore order
  df <- df[order(df$x_pos, df$mutations), ]

  attr(df, "type_order") <- type_order
  attr(df, "medians") <- medians

  return(df)
}


#' Create interactive hamburger plot with plotly
#'
#' Creates a plotly figure with click events for sample selection.
#'
#' @param df Data frame from get_hamburger_data().
#' @param signature_name Character: signature name for title.
#' @param bg_colors Character vector of length 2 for alternating backgrounds.
#' @return A plotly object.
create_interactive_hamburger <- function(
  df,
  signature_name,
  bg_colors = c("#EDF8B1", "#2C7FB8")
) {
  if (is.null(df) || nrow(df) == 0) {
    return(NULL)
  }

  type_order <- attr(df, "type_order")
  medians <- attr(df, "medians")
  n_types <- length(type_order)

  # Create alternating background shapes
  bg_shapes <- lapply(seq_len(n_types), function(i) {
    list(
      type = "rect",
      x0 = i - 0.5,
      x1 = i + 0.5,
      y0 = 0,
      y1 = 1,
      xref = "x",
      yref = "paper",
      fillcolor = bg_colors[(i - 1) %% 2 + 1],
      opacity = 0.2,
      line = list(width = 0),
      layer = "below"
    )
  })

  # Create median line shapes
  median_shapes <- lapply(seq_len(nrow(medians)), function(i) {
    x_pos <- which(type_order == medians$cancer_type[i])
    list(
      type = "line",
      x0 = x_pos - 0.4,
      x1 = x_pos + 0.4,
      y0 = medians$median_mutations[i],
      y1 = medians$median_mutations[i],
      xref = "x",
      yref = "y",
      line = list(color = "red", width = 2, dash = "dash"),
      layer = "above"
    )
  })

  # Combine shapes
  all_shapes <- c(bg_shapes, median_shapes)

  # Create hover text
  df$hover_text <- paste0(
    "<b>Sample:</b> ", df$sample_id, "<br>",
    "<b>Cancer Type:</b> ", df$cancer_type, "<br>",
    "<b>Mutations:</b> ", round(df$mutations, 1)
  )

  # Create the plotly figure
  p <- plotly::plot_ly(
    data = df,
    x = ~x_plot,
    y = ~mutations,
    type = "scatter",
    mode = "markers",
    marker = list(
      size = 6,
      color = "black",
      opacity = 0.6
    ),
    text = ~hover_text,
    hoverinfo = "text",
    customdata = ~sample_id,
    source = "hamburger_plot"
  ) |>
    plotly::layout(
      title = list(
        text = paste0(signature_name, " mutations"),
        x = 0.5
      ),
      xaxis = list(
        title = "",
        tickmode = "array",
        tickvals = seq_len(n_types),
        ticktext = as.character(type_order),
        tickangle = 60,
        side = "top",
        showgrid = FALSE
      ),
      yaxis = list(
        title = "Number of Mutations",
        type = "log",
        showgrid = TRUE
      ),
      shapes = all_shapes,
      hovermode = "closest",
      margin = list(t = 120, b = 50)
    ) |>
    plotly::config(
      displayModeBar = TRUE,
      displaylogo = FALSE
    )

  return(p)
}


#' Get spectrum data for a sample across all ID types
#'
#' Extracts spectrum data for a given sample from 83, 89, and 476 type catalogs.
#'
#' @param sample_id Character: the sample identifier.
#' @param spectra_83 Data frame of ID83 spectra (mutation types as rows).
#' @param spectra_89 Data frame of ID89 spectra.
#' @param spectra_476 Data frame of ID476 spectra.
#' @return List with spectrum data for each ID type.
get_sample_spectra <- function(sample_id, spectra_83, spectra_89, spectra_476) {
  result <- list()

  if (sample_id %in% colnames(spectra_83)) {
    result$id83 <- list(
      mutation_types = rownames(spectra_83),
      counts = as.numeric(spectra_83[, sample_id])
    )
  }

  if (sample_id %in% colnames(spectra_89)) {
    result$id89 <- list(
      mutation_types = rownames(spectra_89),
      counts = as.numeric(spectra_89[, sample_id])
    )
  }

  if (sample_id %in% colnames(spectra_476)) {
    result$id476 <- list(
      mutation_types = rownames(spectra_476),
      counts = as.numeric(spectra_476[, sample_id])
    )
  }

  return(result)
}


#' Define mutation class colors for spectrum plots
#'
#' Returns color mapping consistent with mSigPlot style.
#'
#' @param id_type Character: one of "ID83", "ID89", or "ID476".
#' @return Named list with mutation class colors.
get_mutation_colors <- function(id_type = "ID83") {
  # Colors based on ICAMS/mSigPlot conventions for indel signatures
  # Main categories: DEL (deletions) and INS (insertions)
  # Subcategories by context (C, T, repeat length, etc.)

  colors <- list(
    # 1bp C deletions
    "DEL:C" = "#FDBE6F",
    # 1bp T deletions
    "DEL:T" = "#FC8D59",
    # 2+ bp deletions at repeats
    "DEL:repeats" = "#E34A33",
    # Microhomology deletions
    "DEL:MH" = "#B30000",
    # 1bp C insertions
    "INS:C" = "#78C679",
    # 1bp T insertions
    "INS:T" = "#31A354",
    # 2+ bp insertions at repeats
    "INS:repeats" = "#006837",
    # Complex insertions
    "INS:other" = "#00441B"
  )

  return(colors)
}


#' Classify mutation type into color category
#'
#' @param mut_type Character: mutation type string (e.g., "DEL:C:1:0").
#' @return Character: color category key.
classify_mutation <- function(mut_type) {
  if (grepl("^DEL:C:1:", mut_type)) return("DEL:C")
  if (grepl("^DEL:T:1:", mut_type)) return("DEL:T")
  if (grepl("^DEL:[CT]:[2-9]|^DEL:[CT]:[0-9]{2}", mut_type)) return("DEL:repeats")
  if (grepl("^DEL:.*:M:", mut_type) || grepl("^DEL:MH", mut_type)) return("DEL:MH")
  if (grepl("^INS:C:1:", mut_type)) return("INS:C")
  if (grepl("^INS:T:1:", mut_type)) return("INS:T")
  if (grepl("^INS:[CT]:[2-9]|^INS:[CT]:[0-9]{2}", mut_type)) return("INS:repeats")
  if (grepl("^INS:", mut_type)) return("INS:other")
  if (grepl("^DEL:", mut_type)) return("DEL:repeats")
  return("DEL:C")  # fallback
}


#' Convert spectra to JSON for JavaScript consumption
#'
#' Converts all sample spectra to a JSON object for embedding in HTML.
#'
#' @param spectra_83 Data frame of ID83 spectra.
#' @param spectra_89 Data frame of ID89 spectra.
#' @param spectra_476 Data frame of ID476 spectra.
#' @return JSON string with all spectra data.
spectra_to_json <- function(spectra_83, spectra_89, spectra_476) {
  samples <- unique(c(
    colnames(spectra_83),
    colnames(spectra_89),
    colnames(spectra_476)
  ))

  # Get mutation colors for classification
  colors <- get_mutation_colors()

  # Create nested structure
  all_data <- list(
    id83 = list(
      mutation_types = rownames(spectra_83),
      colors = sapply(rownames(spectra_83), function(mt) {
        cat_key <- classify_mutation(mt)
        colors[[cat_key]]
      }),
      samples = as.list(as.data.frame(spectra_83))
    ),
    id89 = list(
      mutation_types = rownames(spectra_89),
      colors = sapply(rownames(spectra_89), function(mt) {
        cat_key <- classify_mutation(mt)
        colors[[cat_key]]
      }),
      samples = as.list(as.data.frame(spectra_89))
    ),
    id476 = list(
      mutation_types = rownames(spectra_476),
      colors = sapply(rownames(spectra_476), function(mt) {
        cat_key <- classify_mutation(mt)
        colors[[cat_key]]
      }),
      samples = as.list(as.data.frame(spectra_476))
    )
  )

  return(jsonlite::toJSON(all_data, auto_unbox = TRUE))
}
