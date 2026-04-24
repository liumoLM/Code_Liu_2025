# make_spectra_and_plot.R — one-off
#
# Read every indel VCF in this folder, build ID83, ID89, and ID476
# spectra catalogs with ICAMS, and write one multi-page PDF of mSigPlot
# plots next to the VCFs. Each page holds one sample's three panels using
# the fig1.R layout (see ../one_page_83_89_476.R).

library(ICAMS)
library(mSigPlot)
library(ggplot2)

# Locate the directory that contains this script so siblings resolve whether
# sourced interactively, via Rscript, or via source(chdir=).
this_script <- tryCatch(
  {
    ca <- commandArgs(trailingOnly = FALSE)
    m <- regmatches(ca, regexpr("(?<=--file=).+", ca, perl = TRUE))
    if (length(m)) {
      m[1]
    } else {
      ofile <- NULL
      for (i in seq_len(sys.nframe())) {
        f <- sys.frame(i)
        if (!is.null(f$ofile) && is.character(f$ofile) && nzchar(f$ofile)) {
          ofile <- f$ofile
          break
        }
      }
      ofile
    }
  },
  error = function(e) NULL
)
script_dir <- if (!is.null(this_script) && nzchar(this_script)) {
  dirname(normalizePath(this_script))
} else {
  getwd()
}

# ---- Read VCFs ------------------------------------------------------------

vcf_files <- sort(list.files(
  script_dir,
  pattern = "\\.vcf\\.gz$",
  full.names = TRUE
))
stopifnot(length(vcf_files) > 0)

# Strip ".vcf.gz" and the "_INDELintersect" suffix for clean sample names.
sample_names <- sub(
  "_INDELintersect$",
  "",
  sub("\\.vcf\\.gz$", "", basename(vcf_files))
)

split <- ICAMS::ReadAndSplitVCFs(
  files = vcf_files,
  variant.caller = "strelka",
  names.of.VCFs = sample_names
)
id_vcfs <- split$ID
names(id_vcfs) <- sample_names

# Drop variants whose POS lies outside the reference chromosome bounds with
# enough margin for AnnotateIDVCF's sequence-context fetch. At least one input
# VCF contains a chr17 variant past the end of GRCh37/hs37d5 chr17.
bsg <- BSgenome::getBSgenome("BSgenome.Hsapiens.1000genomes.hs37d5")
chr_lens <- GenomeInfoDb::seqlengths(bsg)
margin <- 1000L
for (nm in names(id_vcfs)) {
  v <- id_vcfs[[nm]]
  if (nrow(v) == 0) {
    next
  }
  chroms <- as.character(v$CHROM)
  lens <- chr_lens[chroms]
  indel_w <- pmax(nchar(v$REF), nchar(v$ALT))
  ok <- !is.na(lens) &
    v$POS > margin &
    (v$POS + indel_w + margin) < lens
  dropped <- sum(!ok)
  if (dropped > 0L) {
    message(sprintf(
      "%s: dropping %d out-of-bounds variant(s)",
      nm,
      dropped
    ))
  }
  id_vcfs[[nm]] <- v[ok, , drop = FALSE]
}

# ---- Build catalogs -------------------------------------------------------

# Reference is hs37d5 (see VCF headers) — "hg19" is the matching ICAMS label.
cats <- ICAMS::VCFsToIDCatalogs(
  list.of.vcfs = id_vcfs,
  ref.genome = "hg19",
  region = "genome",
  return.annotated.vcfs = TRUE
)

id83 <- cats$catalog
id476 <- cats$catID476
annot <- cats$annotated.vcfs

# ID89 is not produced by VCFsToIDCatalogs; build it from the annotated VCFs
# (which already carry a Koh_89 column).
id89_cols <- lapply(seq_along(annot), function(i) {
  ICAMS::annot_vcf_to_89_catalog(
    annot_vcf = annot[[i]],
    sample_id = names(annot)[i]
  )
})
id89 <- do.call(cbind, id89_cols)

# Keep column ordering consistent across all three catalogs.
present <- intersect(sample_names, colnames(id83))
id83 <- id83[, present, drop = FALSE]
id89 <- id89[, present, drop = FALSE]
id476 <- id476[, present, drop = FALSE]

# ---- Plot (fig1.R layout: one letter-portrait page per sample) ------------

source(file.path(dirname(script_dir), "one_page_83_89_476.R"))

num_peaks <- 4
base_size <- 9.5
page_w <- 8.5 # letter portrait width  (in)
page_h <- 11 # letter portrait height (in)

out_pdf <- file.path(script_dir, "aa_vcfs_spectra.pdf")

cairo_pdf(out_pdf, width = page_w, height = page_h, onefile = TRUE)
on.exit(grDevices::dev.off(), add = TRUE)

first_page <- TRUE
for (s in present) {
  if (!first_page) {
    grid::grid.newpage()
  }
  first_page <- FALSE
  one_page_83_89_476(
    sig_83 = id83[, s, drop = FALSE],
    sig_89 = id89[, s, drop = FALSE],
    sig_476 = id476[, s, drop = FALSE],
    title_83 = paste0(s, " — ID83"),
    title_89 = paste0(s, " — ID89"),
    title_476 = paste0(s, " — ID476"),
    num_peak_labels = num_peaks,
    base_size = base_size,
    page_h = page_h
  )
}

grDevices::dev.off()
on.exit()

message("Wrote: ", out_pdf)
