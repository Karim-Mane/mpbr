#' Drop a of set SNPs from a `SNPdata` object
#'
#' @param snpdata a `SNPdata` object
#' @param snp_to_be_dropped a `data.frame` with 2 columns" "Chrom" and "Pos"
#' @param chrom the chromosome from which loci should be dropped
#' @param start the starting position of the region to be discarded
#' @param end the end position of the region to be discarded
#'
#' @return a `SNPdata` object where the specified SNPs have been removed
#'
#' @examples
#' \dontrun{
#'   snpdata <- drop_snps(
#'     snpdata,
#'     snp_to_be_dropped = NULL,
#'     chrom             = "Pf3D7_05_v3",
#'     start             = 100,
#'     end               = 500
#'   )
#'  }
#' @details when snp_to_be_dropped is not NULL (i.e. the genomic
#'     coordinates of snps to be removed are in a data frame), then the rest of
#'     the arguments can be ignored or set to NULL (chrom = NULL,
#'     start = NULL, end = NULL)
#' @export
#'
drop_snps <- function(snpdata, snp_to_be_dropped = NULL,
                      chrom = NULL, start = NULL, end = NULL) {
  checkmate::assert_class(snpdata, "SNPdata", null.ok = FALSE)
  checkmate::assert_data_frame(snp_to_be_dropped, ncols = 2L, null.ok = TRUE)
  checkmate::assert_vector(chrom, min.len = 1L, any.missing = FALSE,
                           null.ok = TRUE)
  checkmate::assert_vector(start, min.len = 1L, any.missing = FALSE,
                           null.ok = TRUE)
  checkmate::assert_vector(end, min.len = 1L, any.missing = FALSE,
                           null.ok = TRUE)
  
  if (!is.null(snp_to_be_dropped)) {
    snpdata <- drop_snps_using_gc(snpdata, snp_to_be_dropped)
  } else {
    snpdata <- remove_region_from_snpdata(snpdata, chrom, start, end)
  }
  snpdata[["index"]] <- snpdata[["index"]] + 1L
  snpdata
}

#' Drops SNPs from a `SNPdata` object using their genomic coordinates
#'
#' @inheritParams drop_snps
#'
#' @keywords internal
#'
drop_snps_using_gc <- function(snpdata, snp_to_be_dropped) {
  idx <- which(snpdata[["details"]][["Chrom"]] %in%
                 snp_to_be_dropped[["Chrom"]] &
                 snpdata[["details"]][["Pos"]] %in%
                 snp_to_be_dropped[["Pos"]])
  meta   <- snpdata[["meta"]]
  m      <- which(names(snpdata) %in% c("meta", "vcf", "index"))
  fields <- names(snpdata)[-m]
  for (field in fields) {
    tmp              <- snpdata[[field]][-idx, ]
    snpdata[[field]] <- tmp
  }
  Chrom    <- Pos <- NULL # nolint: object_name_linter
  f2c      <- snpdata[["details"]] %>%
    dplyr::select(Chrom, Pos)
  tmp_file <- file.path(dirname(snpdata[["vcf"]]), "tmp.txt")
  data.table::fwrite(f2c, tmp_file, col.names = FALSE, row.names = FALSE,
                     quote = FALSE, sep = "\t", nThread = 4L)
  index    <- ifelse("index" %in% names(snpdata),
                     snpdata[["index"]] + 1L,
                     0L)
  snpdata[["vcf"]] <- remove_snps_from_vcf(snpdata[["vcf"]],
                                           "tmp.txt",
                                           path  = dirname(snpdata[["vcf"]]),
                                           index = index)
  system(sprintf("rm -f %s", tmp_file))
  snpdata[["index"]] <- index
  snpdata[["meta"]]  <- meta
  snpdata
}

#' Remove SNPs with a specified region of the genome.
#'
#' @inheritParams drop_snps
#'
#' @keywords internal
#'
remove_region_from_snpdata <- function(snpdata, chrom, start, end) {
  stopifnot("Please provide the 'chrom', 'start', 'end' of the region to drop" =
              !is.null(chrom), !is.null(start), !is.null(end))
  idx <- which(snpdata[["details"]][["Chrom"]] == chrom &
                 (snpdata[["details"]][["Pos"]] >= start &
                    snpdata[["details"]][["Pos"]] <= end))
  stopifnot("There is no loci overlapping the specified region." =
              length(idx) > 0L)
  meta   <- snpdata[["meta"]]
  m      <- which(names(snpdata) %in% c("meta", "vcf", "index"))
  fields <- names(snpdata)[-m]
  for (field in fields) {
    tmp              <- snpdata[[field]][-idx, ]
    snpdata[[field]] <- tmp
  }
  Chrom    <- Pos <- NULL # nolint: object_name_linter
  f2c      <- snpdata[["details"]] %>%
    dplyr::select(Chrom, Pos)  # nolint: object_name_linter
  tmp_file <- file.path(dirname(snpdata[["vcf"]]), "tmp.txt")
  data.table::fwrite(f2c, tmp_file, col.names = FALSE,
                     row.names = FALSE, quote = FALSE,
                     sep = "\t", nThread = 4L)
  index    <- ifelse("index" %in% names(snpdata),
                     snpdata[["index"]] + 1L,
                     0L)
  snpdata[["vcf"]] <- remove_snps_from_vcf(snpdata[["vcf"]],
                                           "tmp.txt",
                                           path  = dirname(snpdata[["vcf"]]),
                                           index = index)
  system(sprintf("rm -f %s", tmp_file))
  snpdata[["index"]] <- index
  snpdata[["meta"]]  <- meta
  snpdata            <- validate_snpdata(snpdata)
  snpdata
}