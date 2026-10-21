#' Remove a set of SNPs from a <SNPdata> object
#'
#' @param snpdata The input <SNPdata> object
#' @param coordinates A 2 columns data frame with the genomic coordinates of
#'    SNPs to filter out
#' @param chrom A vector of chromosome names containing the SNPs to be dropped.
#'    This is only used for filtering SNPs within specific regions of the
#'    genome.
#' @param start A numeric vector of region's start positions of the region
#'    covering the SNPs to be discarded.
#' @param end A numeric vector of region's end positions of the region covering
#'    the SNPs to be discarded.
#'
#' @returns A <SNPdata> object where the specified SNPs, regions, or chromosomes
#'    have been removed.
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
remove_snps <- function(snpdata, coordinates = NULL,
                        chrom = NULL, start = NULL, end = NULL) {
  checkmate::assert_class(snpdata, "SNPdata", null.ok = FALSE)
  checkmate::assert_data_frame(coordinates, ncols = 2L, null.ok = TRUE)
  checkmate::assert_vector(chrom, min.len = 1L, any.missing = FALSE,
                           null.ok = TRUE)
  checkmate::assert_vector(start, min.len = 1L, any.missing = FALSE,
                           null.ok = TRUE)
  checkmate::assert_vector(end, min.len = 1L, any.missing = FALSE,
                           null.ok = TRUE)
  # filter SNPs based on provided genomic coordinates
  if (!is.null(coordinates)) {
    snpdata <- drop_snps_using_gc(snpdata, coordinates)
    return(snpdata)
  }
  # filter out an entire chromosome
  if (!is.null(chrom) && all(is.null(c(start, end)))) {
    snpdata <- remove_chrom_from_snpdata(snpdata, chrom)
  } else {
    snpdata <- remove_region_from_snpdata(snpdata, chrom, start, end)
  }
  
  return(snpdata)
}

#' Drops SNPs belonging to the specified chromosome from a <SNPdata> object
#'
#' @inheritParams remove_snps snpdata chrom
#'
#' @keywords internal
#'
remove_chrom_from_snpdata <- function(snpdata, chrom) {
  if (!all(chrom %in% snpdata[["details"]][["Chrom"]])) {
    cli::cli_abort(c(
      i = "Specified {qty(chrom)} chromosome name{?s} not found in the input \\\
      {.cls SNPdata} object.",
      x = "Chromosome name must be one of: \\\
      {.val {toString(unique(snpdata[[\"details\"]][[\"Chrom\"]]))}}."
    ))
  }
  idx <- which(snpdata[["details"]][["Chrom"]] %in% chrom)
  meta <- snpdata[["meta"]]
  m <- which(names(snpdata) %in% c("meta", "vcf"))
  fields <- names(snpdata)[-m]
  for (field in fields) {
    snpdata[[field]] <- snpdata[[field]][-idx, ]
  }
  snpdata[["meta"]] <- meta
  
  # re-calculating the percent of missing SNPs for every sample from the 'GT'
  # matrix
  snpdata[["meta"]][["percentage_missing_sites"]] <-
    colSums(is.na(snpdata[["GT"]])) / nrow(snpdata[["GT"]])
  return(snpdata)
}

#' Drops SNPs from a <SNPdata> object using their genomic coordinates
#'
#' @inheritParams remove_snps snpdata coordinates
#'
#' @keywords internal
#'
drop_snps_using_gc <- function(snpdata, coordinates) {
  idx <- which(
    snpdata[["details"]][["Chrom"]] %in% coordinates[["Chrom"]] &
      snpdata[["details"]][["Pos"]] %in% coordinates[["Pos"]]
  )
  meta <- snpdata[["meta"]]
  m <- which(names(snpdata) %in% c("meta", "vcf"))
  fields <- names(snpdata)[-m]
  for (field in fields) {
    snpdata[[field]] <- snpdata[[field]][-idx, ]
  }
  snpdata[["meta"]] <- meta
  
  # re-calculating the percent of missing SNPs for every sample from the 'GT'
  # matrix
  snpdata[["meta"]][["percentage_missing_sites"]] <-
    colSums(is.na(snpdata[["GT"]])) / nrow(snpdata[["GT"]])
  return(snpdata)
}

#' Remove SNPs within a specified region of the genome.
#'
#' @inheritParams remove_snps snpdata chrom start end
#'
#' @keywords internal
#'
remove_region_from_snpdata <- function(snpdata, chrom, start, end) {
  idx <- which(
    snpdata[["details"]][["Chrom"]] == chrom &
      (snpdata[["details"]][["Pos"]] >= start &
         snpdata[["details"]][["Pos"]] <= end)
  )
  if (!length(idx)) {
    cli::cli_abort(
      "There is no SNP overlapping the specified region."
    )
  }
  meta <- snpdata[["meta"]]
  m <- which(names(snpdata) %in% c("meta", "vcf"))
  fields <- names(snpdata)[-m]
  for (field in fields) {
    snpdata[[field]] <- snpdata[[field]][-idx, ]
  }
  snpdata[["meta"]] <- meta
  
  # re-calculating the percent of missing SNPs for every sample from the 'GT'
  # matrix
  snpdata[["meta"]][["percentage_missing_sites"]] <-
    colSums(is.na(snpdata[["GT"]])) / nrow(snpdata[["GT"]])
  
  return(snpdata)
}