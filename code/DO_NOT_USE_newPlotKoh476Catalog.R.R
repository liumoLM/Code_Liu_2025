library(ggplot2)
library(reshape2)
library(dplyr)
library(ggrepel)

#' Plot indel profile in a 476-channel bar plot for single sample
#' @param Koh476.catalog A indel catalogue of a single sample (476 values)
#' @param text_size Size of text
#' @param plot_title Title of the plot
#' @param num_x_labels Number of top peaks to label on x-axis (NULL = all labels)
#' @return A 476-channel indel profile plot
#' @export
newPlotKoh476Catalog <- function(
  Koh476.catalog,
  text_size = 3,
  plot_title = "test",
  num_x_labels = 10,
  label_size = 2
) {
  # Load Koh476_indeltype if not already in environment
  if (!exists("Koh476_indeltype")) {
    source("code/Koh89_Koh476_Plotting_Functions.R")
  }

  # Ensure catalog is a numeric vector
  if (is.data.frame(Koh476.catalog) || is.matrix(Koh476.catalog)) {
    Koh476.catalog <- as.numeric(Koh476.catalog[, 1])
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

  # Create label data for top peaks only
  if (!is.null(num_x_labels) && num_x_labels > 0) {
    top_indices <- order(muts_basis_melt$freq, decreasing = TRUE)[
      1:num_x_labels
    ]
    label_data <- muts_basis_melt[top_indices, ]
  } else {
    label_data <- muts_basis_melt[0, ] # Empty data frame, no labels
  }

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
  blocks$ymin <- max(muts_basis_melt$freq) * 1.25
  blocks$ymax <- max(muts_basis_melt$freq) * 1.37
  blocks$labels <- c(
    "1bp C",
    "1bp T",
    "1bp C",
    "1bp T",
    ">=2bp",
    ">=2bp",
    "Mh",
    "X"
  )
  blocks$cl <- c(
    "black",
    "black",
    "black",
    "black",
    "white",
    "white",
    "white",
    "white"
  )

  indel_mypalette_fill3 <- c("#000000", "#888888", "#DDDDDD")
  entry3 <- table(Koh476_indeltype$Indel3)
  order_entry3 <- c("Del1", "Ins1", "Del2", "Ins2", "DelMH", "Complex")
  entry3 <- entry3[order_entry3]

  blocks3 <- data.frame(
    Type = unique(Koh476_indeltype$Indel3),
    fill = indel_mypalette_fill3,
    xmin = c(0, cumsum(entry3)[-length(entry3)]) + 0.5,
    xmax = cumsum(entry3) + 0.5
  )
  blocks3$ymin <- max(muts_basis_melt$freq) * 1.37
  blocks3$ymax <- max(muts_basis_melt$freq) * 1.49
  blocks3$labels <- c("Deletion", "Insertion", "Del", "Ins", "DelMH", "X")
  blocks3$cl <- c("black", "black", "white", "white", "white", "white")

  indel_mypalette_fill_all <- c(
    "Del1" = "#fe9f38",
    "Ins1" = "#73bf5d",
    "Del2" = "#f14432",
    "Ins2" = "#4a98c9",
    "DelMH" = "#61409b",
    "X" = "black",
    "Del(2,):M(1,)" = "#61409b",
    "Del(2,):R(1,9)" = "#f14432",
    "Del(C)" = "#fdbe6f",
    "Del(T)" = "#ff8001",
    "Ins(2,):R(0,9)" = "#4a98c9",
    "Ins(C)" = "#b0dd8b",
    "Ins(T)" = "#36a12e",
    "Complex" = "black"
  )

  p <- ggplot2::ggplot(
    data = muts_basis_melt,
    ggplot2::aes(x = x_pos, y = freq, fill = Indel)
  ) +
    ggplot2::geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    ggplot2::xlab("Indel Types") +
    ggplot2::ylab("Count") +
    ggplot2::scale_x_continuous(
      breaks = seq_along(indel_positions),
      labels = NULL,
      limits = c(0.5, length(indel_positions) + 0.5)
    ) +
    ggplot2::ggtitle(plot_title) +
    ggplot2::scale_fill_manual(values = indel_mypalette_fill_all) +
    ggplot2::coord_cartesian(ylim = c(0, max(blocks3$ymax)), clip = "off") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 10, colour = "black"),
      legend.position = "none",
      axis.title.x = ggplot2::element_text(size = 15),
      axis.title.y = ggplot2::element_text(size = 15),
      plot.margin = margin(t = 10, r = 10, b = 80, l = 10)
    ) +
    ggplot2::geom_rect(
      data = blocks,
      ggplot2::aes(
        xmin = xmin,
        ymin = ymin,
        xmax = xmax,
        ymax = ymax,
        fill = Type
      ),
      color = NA,
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
    ggplot2::geom_rect(
      data = blocks3,
      ggplot2::aes(
        xmin = xmin,
        ymin = ymin,
        xmax = xmax,
        ymax = ymax,
        fill = Type
      ),
      color = NA,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = blocks3,
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
    ggplot2::scale_colour_manual(
      values = c("black" = "black", "white" = "white")
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
      max.overlaps = 20,
      min.segment.length = 0,
      inherit.aes = FALSE
    )

  return(p)
}
