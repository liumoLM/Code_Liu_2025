# read_all_sigs

data_dir = "../Manuscript_data"

library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(mSigPlot)

sigfiles = dir(data_dir, pattern = "signatures")

sigs = lapply(
  sigfiles,
  \(basename) {
    read.delim(
      file.path(data_dir, basename),
      sep = '\t',
      check.names = FALSE,
      row.names = 1
    )
  }
)

names(sigs) = c("cosmic", "jin", "koh", "t476", "t83", "t89")
sigcosmic = sigs[["cosmic"]]
sigjin = sigs[["jin"]]
sigkoh = sigs[["koh"]]
sig476 = sigs[["t476"]]
sig83 = sigs[["t83"]]
sig89 = sigs[["t89"]]

specfiles = dir(data_dir, pattern = "spectra")
spectra = lapply(
  specfiles,
  \(basename) {
    read.delim(
      file.path(data_dir, basename),
      sep = '\t',
      check.names = FALSE,
      row.names = 1
    )
  }
)
names(spectra) = c("s476", "s83", "s89")
spec476 = spectra[["s476"]]
spec83 = spectra[["s83"]]
spec89 = spectra[["s89"]]

# all83 = cbind(type83_our_sigs, type83_cosmic_sigs, type83_jin_sigs)

# all89 = cbind(type89_our_sigs, type89_koh_sigs)

if (FALSE) {
  us_v_spectra476 = best_matches(
    plotit = \(zzz, title) plot_476(zzz, plot_title = title),
    out_dir = "476_us_v_spectra",
    file.path(data_dir, "Liu_et_al_final_476_type_signatures.tsv"),
    file.path(data_dir, "Liu_et_al_476_type_spectra.tsv")
  )
  write.table(
    us_v_spectra476,
    "us_v_spectra476.tsv",
    sep = '\t',
    row.names = FALSE
  )
}

linktable = read.delim(file.path(data_dir, "89type_to_83type_connection.tsv"))
colnames(linktable) = c("type89", "exemplar", "type83", "type476", "tmp")

us_v_all83 = read.csv("us_vs_all_83.csv") |>
  select(-X, -max_cosine_id, -max_cosine, -euclidean.x, -euclidean.y)

us_v_koh89 = read.csv("us_v_koh_89.csv") |>
  select(-euclidean) |>
  rename(KohID = ID, cosine_v_koh = cosine)

us_v_spectra476 = read.delim("us_v_spectra476.tsv", sep = '\t') |>
  dplyr::mutate(exemplar476 = gsub("..", "::", ID, fixed = TRUE)) |>
  select(-euclidean, -ID) |>
  rename(cosine_v_exemplar476 = cosine) |>
  mutate(signature = paste0(signature, "_476"))
# Warning: had to put the greek alpha and beta back in

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


#### Do Joins

j1 = full_join(linktable, us_v_all83, by = join_by(type83 == signature))

j2 = full_join(j1, us_v_koh89, by = join_by(type89 == signature))

j3 = full_join(j2, us_v_spectra83, join_by(type83 == signature))

j4 = full_join(j3, us_v_spectra89, join_by(type89 == signature))

j5 = full_join(j4, us_v_spectra476, join_by(type476 == signature))

jall = full_join(j5, koh_v_spectra89, join_by(KohID == signature))

View(jall)

## To do, best_matches476; similar to best_matches89

## To do, plot everything

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

  sigs476 = sigs[["t476"]]

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

  px = plot1(sig83, row$type83)
  py = plot1(
    sigcosmic,
    row$cosmicID,
    paste(row$cosmicID, "cosine vs ours = ", row$cosine_v_cosmic)
  )
  pz = plot1(
    sigjin,
    row$jinID,
    paste(row$jinID, "cosine vs ours = ", row$cosine_v_jin)
  )
  grid.arrange(px, py, pz, nrow = 4, ncol = 1)

  p0 = plot1(
    sig89,
    row$type89,
    title = paste("Extracted 89-type signature", row$type89)
  )
  # toplot = sig89[, row$type89, drop = FALSE]
  # p0 = plot_89(toplot, plot_title = paste("Extracted signature", row$type89))

  toplot = spec89[, row$exemplar, drop = FALSE]
  p1 = plot_89(toplot, plot_title = paste("Original exemplar", row$exemplar))

  toplot = spec89[, row$exemplar89, drop = FALSE]
  p2 = plot_89(
    toplot,
    plot_title = paste(
      "Current exemplar",
      row$exemplar89,
      "cosine to sig =",
      row$cosine_v_exemplar89
    )
  )

  toplot = sigkoh[, row$KohID, drop = FALSE]
  p3 = plot_89(
    toplot,
    plot_title = paste(
      "Koh ID",
      row$KohID,
      paste("cosine to extracted =", row$cosine_v_koh)
    )
  )

  toplot = spec89[, row$exemplar_koh, drop = FALSE]
  p4 = plot_89(
    toplot,
    plot_title = paste(
      row$exemplar_koh,
      "cosine to Koh =",
      row$cosine_koh_v_exemplar_koh
    )
  )

  grid.arrange(p0, p1, p2, p3, p4, nrow = 5, ncol = 1)

  toplot = spec476[, row$exemplar89, drop = FALSE]
  p5 = plot_476(
    toplot,
    simplify_labels = FALSE,
    plot_title = paste(row$exemplar89, "(our type-89 exemplar))")
  )

  toplot = spec476[, row$exemplar_koh, drop = FALSE]
  p6 = plot_476(
    toplot,
    simplify_labels = FALSE,
    plot_title = paste(row$exemplar_koh, "(Koh type-89 exemplar))")
  )

  toplot = spec476[, row$exemplar, drop = FALSE]
  p7 = plot_476(
    toplot,
    simplify_labels = FALSE,
    plot_title = paste(row$exemplar, "(our original type-89 exemplar)")
  )

  sigxx = row$type476
  if (sigxx != "") {
    newid = sub("_476", "", sigxx)
    if (!newid %in% colnames(sigs476)) {
      message(newid, " not found in type 476 signatures, skipping")
      return()
    }
    toplot = sigs476[, newid, drop = FALSE]
    p8 = plot_476(
      toplot,
      simplify_labels = FALSE,
      plot_title = paste("Extracted 476-type sig", sigxx)
    )

    if (is.na(row$exemplar476)) {
      message("No type-476 exemplar for", uniqueid, "skipping")
      return()
    }
    toplot = spec476[, row$exemplar476, drop = FALSE]
    p9 = plot_476(
      toplot,
      simplify_labels = FALSE,
      plot_title = paste(
        row$exemplar476,
        "cosine to our 476 extracted = ",
        row$cosine_v_exemplar476
      )
    )

    grid.arrange(p5, p6, p7, nrow = 3, ncol = 1)
    grid.arrange(p8, p9, nrow = 3, ncol = 1)
  } else {
    grid.arrange(p5, p6, nrow = 2, ncol = 1)
  }
  # Put more code here
  dev.off()
}


plot_row(jall[18, ])

for (i in c(13, 14, 17)) {
  plot_row(jall[i, ])
}

foo = lapply(1:nrow(zz), \(row) plot_row(zz[row, ]))

prop_ins_t_476 = function(catalog) {
  # browser()
  rows = grep("Ins(T)", rownames(catalog), fixed = TRUE)
  allsums = colSums(catalog)
  tsums = colSums(catalog[rows, ])
  tsums / allsums
}

prop_ins_t_83 = function(catalog) {
  rows = grep("INS:C:1", rownames(catalog), fixed = TRUE)
  allsums = colSums(catalog)
  tsums = colSums(catalog[rows, ])
  tsums / allsums
}

prop_ins_t_476(sig89[, rep("InsDel9", 2)])

prop_ins_t_476(sig476[, rep("InsDel9", 2)])

prop_ins_t_476(sigkoh[, rep("InD9b", 2)])
prop_ins_t_476(sigkoh[, rep("InD9a", 2)])

prop_ins_t_83(sigcosmic[, rep("ID9", 2)])
prop_ins_t_83(sig83[, rep("C_ID9", 2)])
prop_ins_t_83(sigjin[, rep("jinID9", 2)])
