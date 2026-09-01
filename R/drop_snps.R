#' Remove a set of SNPs from a <SNPdata> object
#'
#' @param snpdata The input <SNPdata> object
#' @param snps A data frame with 2 columns named as `Chrom` and `Pos`
#' @param chrom A vector of chromosome names containing the SNPs to be dropped.
#'    This is only used for filtering SNPs within specific regions of the
#'    genome.
#' @param start A numeric vector of region's start positions of the region
#'    covering the SNPs to be discarded.
#' @param end A numeric vector of region's end positions of the region covering
#'    the SNPs to be discarded.
#'
#' @returns A <SNPdata> object where the specified SNPs have been removed.
#'
#' @details When the value for the `snps` argument is not NULL, then the rest of
#'     the arguments will be ignored.
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
#'  # remove SNPs
#'  idx_snps <- sample(1:nrow(snpdata[["details"]]), size = 5, replace = FALSE)
#'  snpdata <- remove_snps(
#'    snpdata = snpdata,
#'    snps = snpdata[["details"]][idx_snps, c("Chrom", "Pos")]
#'  )
#' }
remove_snps <- function(snpdata, snps = NULL,
                        chrom = NULL, start = NULL, end = NULL) {
  checkmate::assert_class(snpdata, "SNPdata", null.ok = FALSE)
  checkmate::assert_data_frame(snps, ncols = 2L, null.ok = TRUE)
  checkmate::assert_vector(chrom, min.len = 1L, any.missing = FALSE,
                           null.ok = TRUE)
  checkmate::assert_vector(start, min.len = 1L, any.missing = FALSE,
                           null.ok = TRUE)
  checkmate::assert_vector(end, min.len = 1L, any.missing = FALSE,
                           null.ok = TRUE)
  
  if (!is.null(snps)) {
    snpdata <- drop_snps_using_gc(snpdata, snps)
  } else {
    snpdata <- remove_region_from_snpdata(snpdata, chrom, start, end)
  }
  
  # re-calculating the percent of missing SNPs for every sample from the 'GT'
  # matrix
  snpdata[["meta"]][["percentage_missing_sites"]] <-
    colSums(is.na(snpdata[["GT"]])) / nrow(snpdata[["GT"]])
  
  return(snpdata)
}

#' Drops SNPs from a `SNPdata` object using their genomic coordinates
#'
#' @inheritParams remove_snps
#'
#' @keywords internal
#'
drop_snps_using_gc <- function(snpdata, snp_to_be_dropped) {
  idx <- which(
    snpdata[["details"]][["Chrom"]] %in% snp_to_be_dropped[["Chrom"]] &
      snpdata[["details"]][["Pos"]] %in% snp_to_be_dropped[["Pos"]]
  )
  meta   <- snpdata[["meta"]]
  m      <- which(names(snpdata) %in% c("meta", "vcf"))
  fields <- names(snpdata)[-m]
  for (field in fields) {
    snpdata[[field]] <- snpdata[[field]][-idx, ]
  }
  snpdata[["meta"]]  <- meta
  snpdata
}

#' Remove SNPs with a specified region of the genome.
#'
#' @inheritParams remove_snps snpdata chrom start end
#'
#' @keywords internal
#'
remove_region_from_snpdata <- function(snpdata, chrom, start, end) {
  stopifnot("Please provide the 'chrom', 'start', 'end' of the region to drop" =
              !is.null(chrom), !is.null(start), !is.null(end))
  idx <- which(
    snpdata[["details"]][["Chrom"]] == chrom &
      (snpdata[["details"]][["Pos"]] >= start &
         snpdata[["details"]][["Pos"]] <= end)
  )
  stopifnot("There is no loci overlapping the specified region." =
              length(idx) > 0L)
  meta   <- snpdata[["meta"]]
  m      <- which(names(snpdata) %in% c("meta", "vcf"))
  fields <- names(snpdata)[-m]
  for (field in fields) {
    snpdata[[field]] <- snpdata[[field]][-idx, ]
  }
  snpdata[["meta"]]  <- meta
  
  snpdata
}