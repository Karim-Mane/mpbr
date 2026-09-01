#' Calculate minor allele frequency (MAF) at every SNP
#'
#' @param snpdata An object of class <SNPdata>
#' @param genotype A string with the name of the matrix to be used. Default is
#'    "GT". The other possible values are "Phased", "Imputed", "Phased_imputed".
#' @param include_het A Boolean that specifies whether to account for the
#'    heterozygous allele or not. This can only be activated when
#'    `genotype = "GT"` or `genotype = "Imputed"`.
#' @param name A character with the name of the new column that will created to
#'    store the MAF values. Default is `MAF`.
#'
#' @return The input <SNPdata> object with following 2 additional columns in the
#'    **details** table:
#' \enumerate{
#'   \item MAF: minor allele frequency at every each SNPs
#'   \item MAF_allele: the code for the minor allele. Possible values are:
#'      **1**: the alternate allele is the minor allele,
#'      **0**: the reference allele is the minor allele,
#'      **0/1**: the heterozygous allele is the minor allele,
#'      **0=1**: the reference and alternate alleles have the same frequencies,
#'      **0=1=2**: the three alleles have the same frequencies.
#' }
#'
#' @details If `include_het = FALSE`, the mixed alleles will not be considered
#'    in the MAF calculation.
#'
#' @export
#' @examples
#' \dontrun{
#'  # create the SNPdata object
#'  snpdata <- get_snpdata(
#'    vcf_file = system.file("extdata", "test_data.vcf.gz", package = "mpbr"),
#'    meta_file = system.file("extdata", "sample_metadata.RDS",
#'                             package = "mpbr"),
#'    output_dir = getwd(),
#'    gaf = file.path("data", "PlasmoDB-68_Pfalciparum3D7_Curated_GO.gaf.gz"),
#'    gff = file.path("data", "PlasmoDB-68_Pfalciparum3D7.gff")
#'  )
#'  # calculate MAF
#'  snpdata <- compute_maf(
#'    snpdata = snpdata,
#'    genotype = "GT",
#'    include_het = FALSE,
#'    name = NULL
#'  )
#' }
compute_maf <- function(snpdata, genotype = "GT", include_het = FALSE,
                        name = NULL) {
  checkmate::assert_class(snpdata, "SNPdata", null.ok = FALSE)
  checkmate::assert_logical(include_het, any.missing = FALSE, len = 1L,
                            null.ok = FALSE)
  checkmate::assert_character(genotype, any.missing = FALSE, null.ok = FALSE,
                              len = 1L)
  checkmate::assert_character(name, null.ok = TRUE, len = 1L,
                              any.missing = FALSE)
  genotype <- match.arg(
    genotype,
    choices = c("GT", "Phased", "Imputed", "Phased_imputed")
  )
  if (!genotype %in% c("GT", "Imputed") && include_het) {
    cli::cli_abort(c(
      i = "{.emph include_het} can only be {.val TRUE} for unphased genotype \\\
      matrix ({.val GT} or {.val Imputed}).",
      x = "Chosen genotype matrix is not suitable with the value of the \\\
      {.emph include_het} argument."
    ))
  }
  x <- snpdata[[genotype]]
  ref <- rowSums(x == 0L, na.rm = TRUE)
  alt <- rowSums(x == 1L, na.rm = TRUE)
  het <- rowSums(x == 2L, na.rm = TRUE)
  if (include_het) {
    tmp_mat <- cbind(ref, alt, het)
  } else {
    tmp_mat <- cbind(ref, alt)
  }
  res <- t(apply(tmp_mat, 1L, get_maf))

  if (is.null(name)) {
    name <- "MAF"
  }
  snpdata[["details"]][[name]] <- as.numeric(res[, 1L])
  snpdata[["details"]][["MAF_allele"]] <- as.character(res[, 2L])

  return(snpdata)
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