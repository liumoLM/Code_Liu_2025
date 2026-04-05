#' Read an annotated indel VCF file based on a tumor identifier.
#'
#' @param tumor_id A tumor identifier string, optionally prefixed with
#'   \code{"CancerType::"} which is stripped before lookup. After stripping:
#'   \itemize{
#'     \item If the ID starts with \code{"SP"}, searches for a matching VCF
#'           in the PCAWG directory.
#'     \item Otherwise, searches in the FMH directory.
#'   }
#'
#' @return A \code{data.table} with the VCF contents, or \code{NULL} if no
#'   file is found.
read_annotated_vcf <- function(tumor_id) {
  pcawg_dir <- path.expand("~/MEGA/important_mut_sig_data/pcawg_indel_vcfs")
  h_dir <- path.expand("~/MEGA/important_mut_sig_data/fmh-unfiltered_vcfs/")

  # Strip "CancerType::" prefix if present
  id <- sub("^.*::", "", tumor_id)

  if (grepl("^SP", id)) {
    # PCAWG sample
    pattern <- paste0("*", id, "*annotated.indel.vcf.gz")
    search_dir <- pcawg_dir
  } else {
    # FMH sample
    pattern <- paste0("*", id, ".annotated.indel.vcf.gz")
    search_dir <- h_dir
  }

  hits <- Sys.glob(file.path(search_dir, pattern))
  if (length(hits) == 0) {
    warning("No file found for ", id, " in ", search_dir)
    return(NULL)
  }
  if (length(hits) > 1) {
    warning(
      "Multiple files found for ", id, " in ", search_dir, ":\n",
      paste(hits, collapse = "\n")
    )
    return(NULL)
  }
  data.table::fread(hits)
}
