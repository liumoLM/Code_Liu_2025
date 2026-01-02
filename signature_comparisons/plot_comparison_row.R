# lintr::lint("plot_comparison_row.R", linters = source("../.lintr.R")$value)

#' Compute cosine similarity between two column matrices
#' @param a A column matrix or vector
#' @param b A column matrix or vector
#' @return Cosine similarity (0 to 1)
cosine_sim <- function(a, b) {
  a_vec <- as.numeric(a[, 1])
  b_vec <- as.numeric(b[, 1])

  a_norm <- a_vec / sum(a_vec)
  b_norm <- b_vec / sum(b_vec)

  dist_mat <- rbind(a_norm, b_norm)

  suppressMessages(
    philentropy::distance(dist_mat, method = "cosine", test.na = FALSE)[1]
  )
}

plot_comparison_row = function(row, sigs, spectra) {
  sigcosmic = sigs[["cosmic"]]
  sigjin = sigs[["jin"]]
  sigkoh = sigs[["koh"]]
  sig476 = sigs[["t476"]]
  sig83 = sigs[["t83"]]
  sig89 = sigs[["t89"]]
  spec476 = spectra[["s476"]]
  spec83 = spectra[["s83"]]
  spec89 = spectra[["s89"]]

  rr = function(nn) round(nn, 4)

  if (is.na(row$type89)) {
    return()
  }
  if (row$type89 %in% c("InsDel15", "InsDel16")) {
    return()
  }

  uniqueid = paste0(row$type89, "_", row$type83)
  message(uniqueid)

  plot1 = function(catalog, column, title = column) {
    if (!column %in% colnames(catalog)) {
      stop(column, " not in colnames(catalog)")
    }
    thing = catalog[, column, drop = FALSE]
    if (nrow(catalog) == 89) {
      px = plot_89(thing, plot_title = title)
    } else if (nrow(catalog) == 476) {
      px = plot_476(
        thing,
        simplify_labels = FALSE,
        plot_title = title
      )
    } else if (nrow(catalog) == 83) {
      px = plot_83(thing, plot_title = title)
    } else {
      stop("oops")
    }
    return(px)
  }

  # Collect all plots in a list
  plots83 <- list()
  plots89 = list()
  plots476 = list()

  # Page 1: 83-type signatures
  plots83$px <- plot1(sig83, row$type83)

  plots83$py <- plot1(
    sigcosmic,
    row$cosmicID,
    paste(row$cosmicID, "cosine vs ours = ", rr(row$cosine_v_cosmic))
  )

  orig_exemp_83 = spec83[, row$exemplar, drop = FALSE]
  mysig83 = sig83[, row$type83, drop = FALSE]
  tmp_cos = cosine_sim(orig_exemp_83, mysig83)
  plots83$exemplar = plot1(
    spec83,
    row$exemplar,
    paste(
      row$example,
      "table 1 exemplar in 83-type, cos to sig = ",
      rr(tmp_cos)
    )
  )

  best_exemp_83 = spec83[, row$exemplar89, drop = FALSE]
  mysig83 = sig83[, row$type83, drop = FALSE]
  best_exemplar_cos = cosine_sim(best_exemp_83, mysig83)
  plots83$best_exemplar = plot1(
    spec83,
    row$exemplar89,
    paste(
      row$example,
      "best exemplar in 83-type, cos to sig = ",
      rr(best_exemplar_cos)
    )
  )

  if (row$exemplar != row$exemplar89) {
    message(row$exemplar, " -> ", row$exemplar89)
  }

  if (row$cosine_v_jin > 0.89) {
    plots83$pz <- plot1(
      sigjin,
      row$jinID,
      paste(row$jinID, "cosine vs ours = ", rr(row$cosine_v_jin))
    )
  }

  # Page 2: 89-type signatures and exemplars
  plots89$p0 <- plot1(
    sig89,
    row$type89,
    title = paste(row$type89, "as extracted")
  )

  toplot = spec89[, row$exemplar, drop = FALSE]
  mysig89 = sig89[, row$type89, drop = FALSE]
  tmp_cos = cosine_sim(toplot, mysig89)
  plots89$p1 <- plot_89(
    toplot,
    plot_title = paste(
      "table 1 exemplar",
      row$exemplar,
      "cos to sig = ",
      rr(tmp_cos)
    )
  )

  toplot = spec89[, row$exemplar89, drop = FALSE]
  plots89$p2 <- plot_89(
    toplot,
    plot_title = paste(
      row$exemplar89,
      "best exemplar, cos to sig =",
      rr(row$cosine_v_exemplar89)
    )
  )

  if (row$cosine_v_koh > 0.89) {
    toplot = sigkoh[, row$KohID, drop = FALSE]
    plots89$p3 <- plot_89(
      toplot,
      plot_title = paste(
        "Koh ID",
        row$KohID,
        paste("cosine to extracted =", rr(row$cosine_v_koh))
      )
    )

    toplot = spec89[, row$exemplar_koh, drop = FALSE]
    plots89$p4 <- plot_89(
      toplot,
      plot_title = paste(
        row$exemplar_koh,
        "best Koh exemplar; cos to Koh =",
        rr(row$cosine_koh_v_exemplar_koh)
      )
    )
  }

  # 476-type spectra

  toplot = spec476[, row$exemplar, drop = FALSE]
  plots476$p7 <- plot_476(
    toplot,
    simplify_labels = FALSE,
    plot_title = paste(row$exemplar, "(table 1 type-89 exemplar)")
  )

  toplot = spec476[, row$exemplar89, drop = FALSE]
  plots476$p5 <- plot_476(
    toplot,
    simplify_labels = FALSE,
    plot_title = paste(row$exemplar89, "(our type-89 exemplar)")
  )

  if (row$cosine_v_koh > 0.89) {
    toplot = spec476[, row$exemplar_koh, drop = FALSE]
    plots476$p6 <- plot_476(
      toplot,
      simplify_labels = FALSE,
      plot_title = paste(row$exemplar_koh, "(Koh type-89 exemplar)")
    )
  }

  # Page 4: 476-type signatures (if available)
  sigxx = row$type476
  if (!is.null(sigxx) && sigxx != "" && !is.na(sigxx)) {
    newid = sub("_476", "", sigxx)
    if (newid %in% colnames(sig476)) {
      toplot = sig476[, newid, drop = FALSE]
      plots476$p8 <- plot_476(
        toplot,
        simplify_labels = FALSE,
        plot_title = paste("Extracted 476-type sig", sigxx)
      )

      if (!is.na(row$exemplar476)) {
        toplot = spec476[, row$exemplar476, drop = FALSE]
        plots476$p9 <- plot_476(
          toplot,
          simplify_labels = FALSE,
          plot_title = paste(
            row$exemplar476,
            "our 476 exemplar, cos to our 476 extracted = ",
            rr(row$cosine_v_exemplar476)
          )
        )
      }
    } else {
      message(
        newid,
        " not found in type 476 signatures, skipping 476 sig plots"
      )
    }
  }

  # Write all plots to PDF
  cairo_pdf(
    filename = paste0("plots/", uniqueid, ".pdf"),
    onefile = TRUE,
    height = 10,
    width = 7.5
  )

  # Page(s) for 83-type: arrange in pages of 3
  plots83_list <- unname(plots83)
  n_plots83 <- length(plots83_list)
  plots_per_page <- 3
  n_pages83 <- ceiling(n_plots83 / plots_per_page)

  for (page in seq_len(n_pages83)) {
    start_idx <- (page - 1) * plots_per_page + 1
    end_idx <- min(page * plots_per_page, n_plots83)
    page_plots <- plots83_list[start_idx:end_idx]
    do.call(
      grid.arrange,
      c(
        page_plots,
        list(nrow = max(plots_per_page, length(page_plots)), ncol = 1)
      )
    )
  }

  # Page(s) for 89-type: arrange in pages of 4
  plots89_list <- unname(plots89)
  n_plots89 <- length(plots89_list)
  plots_per_page_89 <- 4
  n_pages89 <- ceiling(n_plots89 / plots_per_page_89)

  for (page in seq_len(n_pages89)) {
    start_idx <- (page - 1) * plots_per_page_89 + 1
    end_idx <- min(page * plots_per_page_89, n_plots89)
    page_plots <- plots89_list[start_idx:end_idx]
    do.call(
      grid.arrange,
      c(
        page_plots,
        list(nrow = max(plots_per_page_89, length(page_plots)), ncol = 1)
      )
    )
  }

  # Page(s) for 476-type spectra: arrange in pages of 3
  plots476_list <- unname(plots476)
  n_plots476 <- length(plots476_list)
  plots_per_page_476 <- 3
  n_pages476 <- ceiling(n_plots476 / plots_per_page_476)

  for (page in seq_len(n_pages476)) {
    start_idx <- (page - 1) * plots_per_page_476 + 1
    end_idx <- min(page * plots_per_page_476, n_plots476)
    page_plots <- plots476_list[start_idx:end_idx]
    do.call(
      grid.arrange,
      c(
        page_plots,
        list(nrow = max(plots_per_page_476, length(page_plots)), ncol = 1)
      )
    )
  }

  dev.off()
  if (best_exemplar_cos < 0.9) {
    message("low relationship to 83-type signature ", rr(best_exemplar_cos))
  }
  return(list(type89_id = row$type89, best_exemplar_cos = best_exemplar_cos))
}
