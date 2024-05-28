#' Calculate minor allele frequency (MAF)
#'
#' Uses the `SNPdata` object to calculate the MAF at every loci
#'
#' @param snpdata a `SNPdata` object
#' @param include_het whether to account for the heterozygous allele or not.
#'    This is only used when `mat_name = "GT"`
#' @param mat_name the name of the matrix to use. default is "GT"
#'
#' @return a `SNPdata` object with 2 additional columns in the **details**
#'    table.
#' \enumerate{
#'   \item MAF: minor allele frequency of each SNPs
#'   \item MAF_allele: 1 if the alternate allele is the minor allele.
#'         0 otherwise
#' }
#'
#' @details if `include_het = FALSE`, the mixed allele will not be considered in
#'    the MAF calculation
#'
#' @examples
#' \dontrun{
#'   snpdata <- compute_maf(
#'    snpdata,
#'    include_het = FALSE,
#'    mat_name    = "GT"
#'  )
#'  }
#' @export
#'
compute_maf <- function(snpdata, include_het = FALSE, mat_name = "GT") {
  checkmate::assert_class(snpdata, "SNPdata", null.ok = FALSE)
  checkmate::assert_logical(include_het, any.missing = FALSE, len = 1L,
                            null.ok = FALSE)
  checkmate::assert_character(mat_name, any.missing = FALSE, null.ok = FALSE,
                              len = 1L)
  x   <- snpdata[[mat_name]]
  ref <- rowSums(x == 0L, na.rm = TRUE)
  alt <- rowSums(x == 1L, na.rm = TRUE)
  het <- rowSums(x == 2L, na.rm = TRUE)
  if (include_het) {
    tmp_mat <- cbind(ref, alt, het)
  } else {
    tmp_mat <- cbind(ref, alt)
  }
  res       <- apply(tmp_mat, 1L, get_maf)
  if ("MAF" %in% names(snpdata[["details"]])) {
    new_maf <- paste0("MAF_", mat_name)
    snpdata[["details"]][[new_maf]] <- as.numeric(res[1L, ])
  } else {
    snpdata[["details"]][["MAF"]]        <- as.numeric(res[1L, ])
    snpdata[["details"]][["MAF_allele"]] <-
      as.factor(as.character(as.numeric(round(res[2L, ]))))
    levels(snpdata[["details"]][["MAF_allele"]]) <-
      dplyr::recode_factor(snpdata[["details"]][["MAF_allele"]], REF = "0",
                           ALT = "1", HET = "2", REF_ALT = "3",
                           REF_ALT_HET = "4")
  }
  
  snpdata
}

#' get the minor allele frequency (MAF) and the corresponding allele
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
      allele <- 0L
    } else if (mat[[1L]] > mat[[2L]]) {
      maf    <- mat[[2L]] / sum(mat[[1L]], mat[[2L]])
      allele <- 1L
    } else {
      maf    <- mat[[2L]] / sum(mat[[1L]], mat[[2L]])
      allele <- 3L
    }
  } else {
    if (mat[[1L]] < mat[[2L]]) {
      minor  <- mat[[1L]]
      allele <- 0L
    } else if (mat[[1L]] > mat[[2L]]) {
      minor  <- mat[[2L]]
      allele <- 1L
    } else {
      minor  <- mat[[2L]]
      allele <- 3L
    }
    
    if (minor < mat[[3L]]) {
      maf    <- minor / sum(mat[[1L]], mat[[2L]], mat[[3L]])
      allele <- 3L
    } else if (minor > mat[[3L]]) {
      maf    <- mat[[3L]] / sum(mat[[1L]], mat[[2L]], mat[[3L]])
      allele <- 2L
    } else {
      maf    <- mat[[3L]] / sum(mat[[1L]], mat[[2L]], mat[[3L]])
      allele <- 4L
    }
  }
  
  return(c(maf, allele))
}