# read_all_sigs

data_dir = "../Manuscript_data"

library(tidyr)
library(dplyr)

sigfiles = dir(data_dir, pattern = "signatures")
sigs = lapply(
  sigfiles,
  \(basename) read.delim(file.path(data_dir, basename), sep = '\t')
)
names(sigs) = c("cosmic", "jin", "koh", "t476", "t83", "t89")


specfiles = dir(data_dir, pattern = "spectra")
spectra = lapply(
  specfiles,
  \(basename) {
    read.delim(file.path(data_dir, basename), sep = '\t', check.names = FALSE)
  }
)
names(spectra) = c("s476", "s83", "s89")

# all83 = cbind(type83_our_sigs, type83_cosmic_sigs, type83_jin_sigs)

# all89 = cbind(type89_our_sigs, type89_koh_sigs)

if (FALSE) {
  us_v_spectra476 = best_matches(
    plotit = \(zzz, title) plot_476(zzz, plot_title = title),
    out_dir = "476_us_v_spectra",
    file.path(data_dir, "Liu_et_al_final_476_type_signatures.tsv"),
    file.path(data_dir, "Liu_et_al_476_type_spectra.tsv")
  )
  write.csv(us_v_spectra476, "us_v_spectra476.csv")
}

linktable = read.delim(
  "linktable.tsv",
  sep = '\t'
)

us_v_all83 = read.csv(
  "us_vs_all_83.csv"
) |>
  select(-X, -max_cosine_id, -max_cosine, -euclidean.x, -euclidean.y)

us_v_koh89 = read.csv(
  "us_v_koh_89.csv"
) |>
  select(-euclidean)

us_v_spectra476 = read.csv(
  "us_v_spectra476.csv"
) |>
  dplyr::mutate(exemplar476 = gsub("..", "::", ID, fixed = TRUE)) |>
  select(-euclidean, -X, -ID) |>
  rename(cosine_v_exemplar476 = cosine) |>
  mutate(signature = paste0(signature, "_476"))

us_v_spectra83 = read.delim("us_v_spectra_83.tsv", sep = '\t') |>
  dplyr::mutate(exemplar83 = gsub("..", "::", ID, fixed = TRUE)) |>
  select(-euclidean, -ID) |>
  rename(cosine_v_exemplar83 = cosine)

us_v_spectra89 = read.delim("us_v_spectra_89.tsv", sep = '\t') |>
  dplyr::mutate(exemplar89 = gsub("..", "::", ID, fixed = TRUE)) |>
  select(-euclidean, -ID) |>
  rename(cosine_v_exemplar89 = cosine)

koh_v_spectra89 = read.delim("koh_v_spectra_89.tsv") |>
  dplyr::mutate(exemplar_koh = gsub("..", "::", ID, fixed = TRUE)) |>
  select(-euclidean, -ID) |>
  rename(cosine_koh_v_exemplar_koh = cosine)


#### Start joins

link_and_83 = dplyr::full_join(
  linktable,
  us_v_all83,
  by = join_by(type83 == signature)
)

link_and_83_and_89 = dplyr::full_join(
  link_and_83,
  us_v_koh89,
  by = join_by(type89 == signature)
) |>
  rename(KohID = ID, cosine_v_koh = cosine)


link_and_83_and_89_and_spectra83 = dplyr::full_join(
  link_and_83_and_89,
  us_v_spectra83,
  join_by(type83 == signature)
)


link_and_83_and_89_and_spectra83_and_spectra89 = dplyr::full_join(
  link_and_83_and_89_and_spectra83,
  us_v_spectra89,
  join_by(type89 == signature)
)

link_and_83_and_89_and_all_type_spectra = dplyr::full_join(
  link_and_83_and_89_and_spectra83_and_spectra89,
  us_v_spectra476,
  join_by(type476 == signature)
)

moremore = full_join(
  link_and_83_and_89_and_all_type_spectra,
  koh_v_spectra89,
  join_by(KohID == signature)
)

View(moremore)

## To do, best_matches476; similar to best_matches89

## To do, plot everything

source("../code/wrap_ICAMS_plot_catalog.R")

plot_row = function(row) {
  # browser()

  if (is.na(row$type89)) {
    return()
  }
  if (row$type89 %in% c("InsDel15", "InsDel16")) {
    return()
  }

  uniqueid = paste0(row$type89, "_", row$type83)
  message(uniqueid)
  cairo_pdf(
    filename = paste0("plots/", uniqueid, ".pdf"),
    onefile = TRUE,
    height = 10,
    width = 7.5
  )
  par(mfrow = c(5, 1))
  # browser()
  s89 = spectra[["s89"]]
  s476 = spectra[["s476"]]
  koh = sigs[["koh"]]
  sigs89 = sigs[["t89"]]
  sigs476 = sigs[["t476"]]

  toplot = sigs89[, row$type89, drop = FALSE]
  p0 = plot_89(toplot, plot_title = paste("Extracted signature", row$type89))

  toplot = s89[, row$exemplar, drop = FALSE]
  p1 = plot_89(toplot, plot_title = paste("Orignal exemplar", row$exemplar))

  toplot = s89[, row$exemplar89, drop = FALSE]
  p2 = plot_89(
    toplot,
    plot_title = paste(
      "Current exemplar",
      row$exemplar,
      "cosine to sig =",
      row$cosine_v_exemplar89
    )
  )

  toplot = koh[, row$KohID, drop = FALSE]
  p3 = plot_89(
    toplot,
    plot_title = paste(
      "Koh ID",
      row$KohID,
      paste("cosine to extracted =", row$cosine_v_koh)
    )
  )

  toplot = s89[, row$exemplar_koh, drop = FALSE]
  p4 = plot_89(
    toplot,
    plot_title = paste(
      row$exemplar_koh,
      "cosine to Koh =",
      row$cosine_koh_v_exemplar_koh
    )
  )

  grid.arrange(p0, p1, p2, p3, p4, nrow = 5, ncol = 1)

  toplot = s476[, row$exemplar89, drop = FALSE]
  p5 = plot_476(
    toplot,
    simplify_labels = FALSE,
    plot_title = paste(row$exemplar89, "(our type-89 exemplar))")
  )

  toplot = s476[, row$exemplar_koh, drop = FALSE]
  p6 = plot_476(
    toplot,
    simplify_labels = FALSE,
    plot_title = paste(row$exemplar_koh, "(Koh type-89 exemplar))")
  )

  sigxx = row$type476
  if (sigxx != "") {
    newid = sub("_476", "", sigxx)
    if (!newid %in% colnames(sigs476)) {
      message(newid, " not found in type 476 signatures, skipping")
      return()
    }
    toplot = sigs476[, newid, drop = FALSE]
    p7 = plot_476(
      toplot,
      simplify_labels = FALSE,
      plot_title = paste("Extracted 476-type sig", sigxx)
    )

    if (is.na(row$exemplar476)) {
      message("No type-476 exemplar for", uniqueid, "skipping")
      return()
    }
    toplot = s476[, row$exemplar476, drop = FALSE]
    p8 = plot_476(
      toplot,
      simplify_labels = FALSE,
      plot_title = paste(
        row$exemplar476,
        "cosine to our 476 extracted = ",
        row$cosine_v_exemplar476
      )
    )

    grid.arrange(p5, p6, p7, nrow = 3, ncol = 1)
    grid.arrange(p8, nrow = 3, ncol = 1)
  } else {
    grid.arrange(p5, p6, nrow = 2, ncol = 1)
  }
  # Put more code here
  dev.off()
}

zz = moremore

plot_row(zz[1, ])

foo = lapply(1:nrow(zz), \(row) plot_row(zz[row, ]))
