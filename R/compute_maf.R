#' Calculate minor allele frequency (MAF) at every loci
#'
#' @param snpdata An object of class `SNPdata`
#' @param mat_name A string with the name of the matrix to be used. Default is
#'    "GT". The other possible values are "Phased", "Imputed", "Phased_imputed".
#' @param include_het A Boolean that specifies whether to account for the
#'    heterozygous allele or not. This can only be activated when
#'    `mat_name = "GT"`.
#'
#' @return The input `SNPdata` object with following 2 additional columns in the
#'    **details** table:
#' \enumerate{
#'   \item MAF: minor allele frequency at every each SNPs
#'   \item MAF_allele: the code for the minor allele. Possible values are:
#'      **1**: the alternate allele is the minor allele,
#'      **0**: the alternate allele is the minor allele,
#'      **0/1**: the heterozygous allele is the minor allele,
#'      **0=1**: the reference and alternate alleles have the same frequencies,
#'      **0=1=2**: the three alleles have the same frequencies.
#' }
#'
#' @details If `include_het = FALSE`, the mixed alleles will not be considered
#'    in the MAF calculation.
#'
#' @examples
#' \dontrun{
#'   snpdata <- compute_maf(
#'     snpdata,
#'     include_het = FALSE,
#'     mat_name    = "GT"
#'   )
#' }
#' @export
#'
compute_maf <- function(snpdata, include_het = FALSE, mat_name = "GT") {
  checkmate::assert_class(snpdata, "SNPdata", null.ok = FALSE)
  checkmate::assert_logical(include_het, any.missing = FALSE, len = 1L,
                            null.ok = FALSE)
  checkmate::assert_character(mat_name, any.missing = FALSE, null.ok = FALSE,
                              len = 1L)
  mat_name <- match.arg(
    mat_name,
    choices = c("GT", "Phased", "Imputed", "Phased_imputed")
  )
  x   <- snpdata[[mat_name]]
  ref <- rowSums(x == 0L, na.rm = TRUE)
  alt <- rowSums(x == 1L, na.rm = TRUE)
  het <- rowSums(x == 2L, na.rm = TRUE)
  if (include_het) {
    tmp_mat <- cbind(ref, alt, het)
  } else {
    tmp_mat <- cbind(ref, alt)
  }
  res       <- t(apply(tmp_mat, 1L, get_maf))
  if ("MAF" %in% names(snpdata[["details"]])) {
    new_maf <- paste0("MAF_", mat_name)
    snpdata[["details"]][[new_maf]] <- as.numeric(res[, 1L])
  } else {
    snpdata[["details"]][["MAF"]]        <- as.numeric(res[, 1L])
    snpdata[["details"]][["MAF_allele"]] <- as.character(res[, 2L])
  }
  
  snpdata
}

#' Get the minor allele frequency (MAF) and the corresponding allele
#'
#' @param mat a matrix with the 2 or 3 columns. Every column should contain the
#'    the count of the allele of interest across all samples.
#'
#' @return a `vector` of the MAF and the corresponding allele
#' @keywords internal
#'
get_maf <- function(mat) {
  if (length(mat) == 2L) {
    if (mat[[1L]] < mat[[2L]]) {
      maf    <- mat[[1L]] / sum(mat[[1L]], mat[[2L]])
      allele <- "0"
    } else if (mat[[1L]] > mat[[2L]]) {
      maf    <- mat[[2L]] / sum(mat[[1L]], mat[[2L]])
      allele <- "1"
    } else {
      maf    <- mat[[2L]] / sum(mat[[1L]], mat[[2L]])
      allele <- "0=1"
    }
  } else {
    if (mat[[1L]] < mat[[2L]]) {
      minor  <- mat[[1L]]
      allele <- "0"
    } else if (mat[[1L]] > mat[[2L]]) {
      minor  <- mat[[2L]]
      allele <- "1"
    } else {
      minor  <- mat[[2L]]
      allele <- "0=1"
    }
    
    if (minor < mat[[3L]]) {
      maf    <- minor / sum(mat[[1L]], mat[[2L]], mat[[3L]])
      # allele <- 3L
    } else if (minor > mat[[3L]]) {
      maf    <- mat[[3L]] / sum(mat[[1L]], mat[[2L]], mat[[3L]])
      allele <- "0/1"
    } else {
      maf    <- mat[[3L]] / sum(mat[[1L]], mat[[2L]], mat[[3L]])
      allele <- "0=1=2"
    }
  }
  
  return(c(maf, allele))
}