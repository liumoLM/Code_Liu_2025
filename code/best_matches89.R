library(gridExtra)

source("code/Koh89_Koh476_Plotting_Functions.R")
source("code/best_matches.R")

#' Run best_matches for ID89 signatures
#' @param out_dir Output directory for PDF files
#' @param ... At least 2 file paths: first is our_sigs, rest are reference sigs
best_matches89 <- function(out_dir, ...) {
  paths <- list(...)
  if (length(paths) < 2) {
    stop("At least 2 file paths required: our_sigs and at least one reference")
  }

  # Define plotit function using PlotKoh89Catalog directly
  plotit <- function(vec, title) {
    PlotKoh89Catalog(vec, plot_title = title)
  }

  best_matches(
    plotit,
    out_dir,
    ...
  )
}


us_v_koh_89 = best_matches89(
  "89_us_v_koh",
  "data/type89_our_sigs.tsv",
  "data/type89_koh_sigs.tsv"
)

koh_v_us_89 =
  best_matches89(
    "89_koh_v_us",
    "data/type89_koh_sigs.tsv",
    "data/type89_our_sigs.tsv"
  )

IndelType = c(
  # Single base deletions - C
  "[Del(C):R1]A",
  "[Del(C):R1]T",
  "[Del(C):R2]A",
  "[Del(C):R2]T",
  "[Del(C):R3]A",
  "[Del(C):R3]T",
  "[Del(C):R(4,5)]A",
  "[Del(C):R(4,5)]T",
  "[Del(C):R(1,5)]G",
  "Del(C):R(6,9)",
  # Single base deletions - T
  "A[Del(T):R(1,4)]A",
  "A[Del(T):R(1,4)]C",
  "A[Del(T):R(1,4)]G",
  "C[Del(T):R(1,4)]A",
  "C[Del(T):R(1,4)]C",
  "C[Del(T):R(1,4)]G",
  "G[Del(T):R(1,4)]A",
  "G[Del(T):R(1,4)]C",
  "G[Del(T):R(1,4)]G",
  "A[Del(T):R(5,7)]A",
  "A[Del(T):R(5,7)]C",
  "A[Del(T):R(5,7)]G",
  "C[Del(T):R(5,7)]A",
  "C[Del(T):R(5,7)]C",
  "C[Del(T):R(5,7)]G",
  "G[Del(T):R(5,7)]A",
  "G[Del(T):R(5,7)]C",
  "G[Del(T):R(5,7)]G",
  "A[Del(T):R(8,)]A",
  "A[Del(T):R(8,)]C",
  "A[Del(T):R(8,)]G",
  "C[Del(T):R(8,)]A",
  "C[Del(T):R(8,)]C",
  "C[Del(T):R(8,)]G",
  "G[Del(T):R(8,)]A",
  "G[Del(T):R(8,)]C",
  "G[Del(T):R(8,)]G",
  # Single base insertions - C
  "A[Ins(C):R0]A",
  "A[Ins(C):R0]T",
  "Ins(C):R(0,3)",
  "Ins(C):R(4,6)",
  "Ins(C):R(7,)",
  # Single base insertions - T
  "A[Ins(T):R(0,4)]A",
  "A[Ins(T):R(0,4)]C",
  "A[Ins(T):R(0,4)]G",
  "C[Ins(T):R(0,4)]A",
  "C[Ins(T):R(0,4)]C",
  "C[Ins(T):R(0,4)]G",
  "G[Ins(T):R(0,4)]A",
  "G[Ins(T):R(0,4)]C",
  "G[Ins(T):R(0,4)]G",
  "A[Ins(T):R(5,7)]A",
  "A[Ins(T):R(5,7)]C",
  "A[Ins(T):R(5,7)]G",
  "C[Ins(T):R(5,7)]A",
  "C[Ins(T):R(5,7)]C",
  "C[Ins(T):R(5,7)]G",
  "G[Ins(T):R(5,7)]A",
  "G[Ins(T):R(5,7)]C",
  "G[Ins(T):R(5,7)]G",
  "A[Ins(T):R(8,)]A",
  "A[Ins(T):R(8,)]C",
  "A[Ins(T):R(8,)]G",
  "C[Ins(T):R(8,)]A",
  "C[Ins(T):R(8,)]C",
  "C[Ins(T):R(8,)]G",
  "G[Ins(T):R(8,)]A",
  "G[Ins(T):R(8,)]C",
  "G[Ins(T):R(8,)]G",
  # Longer deletions (no MH)
  "Del(2,4):R1",
  "Del(5,):R1",
  "Del(2,8):U(1,2):R(2,4)",
  "Del(2,):U(1,2):R(5,)",
  "Del(3,):U(3,):R2",
  "Del(3,):U(3,):R(3,)",
  # Longer insertions
  "Ins(2,4):R0",
  "Ins(5,):R0",
  "Ins(2,4):R1",
  "Ins(5,):R1",
  "Ins(2,):R(2,4)",
  "Ins(2,):R(5,)",
  # Deletions with MH
  "Del(2,5):M1",
  "Del(3,5):M2",
  "Del(4,5):M(3,4)",
  "Del(6,):M1",
  "Del(6,):M2",
  "Del(6,):M3",
  "Del(6,):M(4,)",
  # Complex
  "Complex"
)

uu = rownames(yy)
aa = gsub(
  "R\\(5,9",
  "R\\(5,",
  gsub("R\\(3,9", "R\\(3,", gsub("R\\(8,9", "R\\(8,", uu))
)
bb = gsub("R\\(7,9", "R\\(7,", aa)
rownames(yy) = bb
dim(yy)
yy2 = yy[IndelType, ]
write.table(yy2, "data/type89_koh_sigs2.tsv", sep = '\t')
zz = read.delim("data/type89_our_sigs.tsv")
rownames(zz) == IndelType
