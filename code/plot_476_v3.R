library(ggplot2)
library(reshape2)
library(dplyr)
library(ggrepel)

#' Plot indel profile in a 476-channel bar plot for single sample
#' with grouped flanking base labels below the x-axis
#' @param Koh476.catalog A indel catalogue of a single sample (476 values)
#' @param text_size Size of text
#' @param plot_title Title of the plot
#' @param num_x_labels Number of top peaks to label per block (NULL = none)
#' @param label_size Size of peak labels
#' @param threshold_denominator Denominator for minimum threshold calculation (default = 7)
#' @return A 476-channel indel profile plot
#' @export
plot_476_v3 <- function(
  Koh476.catalog,
  text_size = 3,
  plot_title = "test",
  num_x_labels = 3,
  label_size = 2,
  threshold_denominator = 7
) {
  # Load Koh476_indeltype if not already in environment
  if (!exists("Koh476_indeltype")) {
    source("code/Koh89_Koh476_Plotting_Functions.R")
  }

  # Ensure catalog is a numeric vector
  if (is.data.frame(Koh476.catalog) || is.matrix(Koh476.catalog)) {
    Koh476.catalog <- as.numeric(Koh476.catalog[, 1])
  }

  # Determine y-axis label based on sum
  ylabel <- if (sum(Koh476.catalog, na.rm = TRUE) < 1.1) {
    "Proportion"
  } else {
    "Count"
  }

  my_vector <- Koh476_indeltype$IndelType
  muts_basis <- data.frame(Sample = Koh476.catalog, IndelType = my_vector)
  muts_basis_melt <- reshape2::melt(muts_basis, "IndelType")
  muts_basis_melt <- merge(
    Koh476_indeltype,
    muts_basis_melt,
    by = "IndelType",
    all.x = TRUE
  )
  muts_basis_melt[is.na(muts_basis_melt)] <- 0
  names(muts_basis_melt) <- c(
    "IndelType",
    "Indel",
    "Indel3",
    "Figlabel",
    "Sample",
    "freq"
  )
  muts_basis_melt$Sample <- as.character(muts_basis_melt$Sample)

  # Add x position for plotting
  muts_basis_melt$x_pos <- match(
    muts_basis_melt$IndelType,
    Koh476_indeltype$IndelType
  )

  indel_mypalette_fill <- c(
    "#000000",
    "#61409b",
    "#f14432",
    "#fdbe6f",
    "#ff8001",
    "#4a98c9",
    "#b0dd8b",
    "#36a12e"
  )

  indel_positions <- Koh476_indeltype$IndelType

  entry <- table(Koh476_indeltype$Indel)
  order_entry <- c(
    "Del(C)",
    "Del(T)",
    "Ins(C)",
    "Ins(T)",
    "Del(2,):R(1,9)",
    "Ins(2,):R(0,9)",
    "Del(2,):M(1,)",
    "Complex"
  )
  entry <- entry[order_entry]

  blocks <- data.frame(
    Type = unique(Koh476_indeltype$Indel),
    fill = indel_mypalette_fill,
    xmin = c(0, cumsum(entry)[-length(entry)]) + 0.5,
    xmax = cumsum(entry) + 0.5
  )
  blocks$ymin <- max(muts_basis_melt$freq) * 1.35
  blocks$ymax <- max(muts_basis_melt$freq) * 1.47
  blocks$labels <- c(
    "Del 1bp C",
    "Del 1bp T",
    "Ins 1bp C",
    "Ins 1bp T",
    "Del ≥2bp",
    "Ins ≥2bp",
    "Mh",
    "X"
  )
  blocks$cl <- c(
    "black",
    "black",
    "black",
    "black",
    "black",
    "black",
    "black",
    "white"
  )

  indel_mypalette_fill_all <- c(
    "Del(2,):M(1,)" = "#61409b",
    "Del(2,):R(1,9)" = "#f14432",
    "Del(C)" = "#fdbe6f",
    "Del(T)" = "#ff8001",
    "Ins(2,):R(0,9)" = "#4a98c9",
    "Ins(C)" = "#b0dd8b",
    "Ins(T)" = "#36a12e",
    "Complex" = "black"
  )

  # Create label data: top num_x_labels per block, excluding peaks < 1/10 max
  max_freq <- max(muts_basis_melt$freq)
  min_threshold <- max_freq / threshold_denominator

  if (!is.null(num_x_labels) && num_x_labels > 0) {
    label_data <- muts_basis_melt |>
      dplyr::filter(freq >= min_threshold) |>
      dplyr::group_by(Indel) |>
      dplyr::slice_max(order_by = freq, n = num_x_labels, with_ties = FALSE) |>
      dplyr::ungroup()
  } else {
    label_data <- muts_basis_melt[0, ] # Empty data frame, no labels
  }

  # Create flanking base annotation blocks
  flanking_patterns <- c(
    "A\\[Del",
    "G\\[Del",
    "T\\[Del",
    "A\\[Ins",
    "G\\[Ins",
    "T\\[Ins",
    "C\\[Del",
    "C\\[Ins"
  )
  flanking_labels <- c(
    "⎢A...",
    "⎢G...",
    "⎢T...",
    "⎢A...",
    "⎢G...",
    "⎢T...",
    "⎢C...",
    "⎢C..."
  )

  # Find runs of each pattern
  flanking_blocks <- data.frame(
    label = character(),
    xmin = numeric(),
    xmax = numeric(),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(flanking_patterns)) {
    matches <- grep(flanking_patterns[i], Koh476_indeltype$IndelType)
    if (length(matches) > 0) {
      # Find contiguous runs
      runs <- split(matches, cumsum(c(1, diff(matches) != 1)))
      for (run in runs) {
        flanking_blocks <- rbind(
          flanking_blocks,
          data.frame(
            label = flanking_labels[i],
            xmin = min(run) - 0.5,
            xmax = max(run) + 0.5,
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }

  # Position flanking blocks below the x-axis
  flanking_blocks$y <- -max(muts_basis_melt$freq) * 0.17

  # Find positions for vertical lines
  pos_a_del_c <- which(Koh476_indeltype$IndelType == "A[Del(C):R1]A")
  pos_g_del_c <- which(Koh476_indeltype$IndelType == "G[Del(C):R1]A")

  p <- ggplot2::ggplot(
    data = muts_basis_melt,
    ggplot2::aes(x = x_pos, y = freq, fill = Indel)
  ) +
    ggplot2::geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    ggplot2::geom_vline(
      xintercept = pos_a_del_c,
      linetype = "dashed",
      color = "gray50"
    ) +
    ggplot2::geom_vline(
      xintercept = pos_g_del_c,
      linetype = "dashed",
      color = "gray50"
    ) +
    ggplot2::xlab("Indel Type") +
    ggplot2::ylab(ylabel) +
    ggplot2::scale_x_continuous(
      breaks = NULL,
      labels = NULL,
      limits = c(0.5, length(indel_positions) + 0.5)
    ) +
    ggplot2::ggtitle(plot_title) +
    ggplot2::scale_fill_manual(values = indel_mypalette_fill_all) +
    ggplot2::coord_cartesian(ylim = c(0, max(blocks$ymax)), clip = "off") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 10, colour = "black"),
      legend.position = "none",
      axis.title.x = ggplot2::element_text(size = 10, margin = margin(t = 10)),
      axis.title.y = ggplot2::element_text(size = 10),
      plot.margin = margin(t = 10, r = 10, b = 80, l = 10)
    ) +
    ggplot2::scale_colour_manual(
      values = c("black" = "black", "white" = "white")
    ) +
    ggplot2::geom_rect(
      data = blocks,
      ggplot2::aes(
        xmin = xmin,
        ymin = ymin,
        xmax = xmax,
        ymax = ymax,
        fill = Type,
        colour = "white"
      ),
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = blocks,
      ggplot2::aes(
        x = (xmax + xmin) / 2,
        y = (ymax + ymin) / 2,
        label = labels,
        colour = cl
      ),
      size = text_size,
      fontface = "bold",
      inherit.aes = FALSE
    ) +
    # Add flanking base labels below x-axis
    ggplot2::geom_text(
      data = flanking_blocks,
      ggplot2::aes(
        x = xmin + 6,
        y = y,
        label = label
      ),
      size = 2,
      inherit.aes = FALSE
    ) +
    ggrepel::geom_text_repel(
      data = label_data,
      ggplot2::aes(x = x_pos, y = freq, label = Figlabel),
      size = label_size,
      nudge_y = max(muts_basis_melt$freq) * 0.1,
      direction = "both",
      segment.color = "gray50",
      segment.size = 0.3,
      arrow = grid::arrow(length = unit(0.02, "npc"), type = "closed"),
      max.overlaps = 50,
      min.segment.length = 0,
      inherit.aes = FALSE
    )

  return(p)
}
