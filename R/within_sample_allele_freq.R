#' Compute the within-sample minor allele frequency (WSMAF)
#'
#' @param vcf_file A character with the path to the input VCF file
#'
#' @returns A numeric matrix with the within-sample minor allele frequency
#'    values across all loci
#' @export
#'
#' @examples
#' \dontrun{
#'   wsmaf <- compute_wsmaf(
#'     vcf_file = system.file("extdata", "test_data.vcf.gz", package = "mpbr")
#'   )
#' }
compute_wsmaf <- function(vcf_file) {
  checkmate::assert_file_exists(vcf_file, access = "r")

  # import the VCF file
  vcf_data <- data.table::fread(
    file = vcf_file,
    skip = "#CHROM",
    header = TRUE,
    sep = "\t",
    nThread = 5
  )

  # extract the allelic depth
  allelic_depth <- mpbr::extract_genotype(
    vcf_data = vcf_data,
    which = "AD",
    nthreads = 10
  )

  # compute the WSMAF
  names(allelic_depth)[1] <- "CHROM"
  tmp <- allelic_depth |>
    dplyr::select(-c(CHROM, POS))
  wsmaf <- apply(tmp, MARGIN = c(1, 2), calculate_af)
  wsmaf <- apply(wsmaf, 2, as.numeric)
  rownames(wsmaf) <- paste0(
    allelic_depth[["CHROM"]], "_", allelic_depth[["POS"]]
  )
  return(t(wsmaf))
}

#' Calculate the within-sample minor allele frequency at a specific locus
#'
#' @param x A character with the allelic depth of a sample at a specific SNP
#'    position
#'
#' @returns A numeric value corresponding to the within-sample minor allele
#'    frequency at a specific locus
#' @keywords internal
calculate_af <- function(x) {
  if (x == ".") {
    return(NA)
  }
  x <- as.numeric(unlist(strsplit(x, ",", fixed = TRUE)))
  x <- x / sum(x)
  return(as.numeric(x[which.min(x)]))
}