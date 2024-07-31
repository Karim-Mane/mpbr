#' Filter a `SNPdata` object
#'
#' SNPs and samples are filtered based on the user-specified thresholds. Both
#' the `meta`, `details` and the different genotype tables will be affected
#' depending on the filtration conditions.
#'
#' @param snpdata An object of class `SNPdata`
#' @param min_qual An integer that represents the The minimum call quality
#'    score. Loci with a quality score below this will be discarded.
#'    Default is 1000.
#' @param max_missing_sites A numeric representing the maximum fraction of
#'    missing sites above which a sample should be dropped. Default is 0.2.
#' @param max_missing_samples A numeric representing the maximum fraction of
#'    missing samples above which a loci should be discarded. Default is 0.2.
#' @param maf_cutoff A numeric representing the minor allele frequency cut-off.
#'    Loci with a `MAF < maf_cutoff` will be discarded.
#'
#' @return a filtered `SNPdata` object
#' @examples
#' \dontrun{
#'  snpdata <- filter(
#'   snpdata,
#'   min_qual            = 1000,
#'   max_missing_sites   = 0.2,
#'   max_missing_samples = 0.2,
#'   maf_cutoff          = 0.01
#'  )
#'  }
#' @export
#'
filter <- function(snpdata,
                   min_qual            = 1000L,
                   max_missing_sites   = 0.2,
                   max_missing_samples = 0.2,
                   maf_cutoff          = 0.01) {
  checkmate::assert_class(snpdata, "SNPdata", null.ok = FALSE)
  checkmate::assert_numeric(min_qual, lower = 1000L, any.missing = FALSE,
                            null.ok = FALSE)
  checkmate::assert_numeric(max_missing_sites, lower = 0L, upper = 1L,
                            finite = TRUE, any.missing = FALSE, null.ok = FALSE,
                            len = 1L)
  checkmate::assert_numeric(max_missing_samples, lower = 0L, upper = 1L,
                            finite = TRUE, any.missing = FALSE, null.ok = FALSE,
                            len = 1L)
  checkmate::assert_numeric(maf_cutoff, lower = 0L, upper = 1L,
                            finite = TRUE, any.missing = FALSE, null.ok = FALSE,
                            len = 1L)
  # filter the details table and the genotype matrices
  x      <- snpdata[["details"]]
  fields <- names(snpdata)[!(names(snpdata) %in% c("meta", "details", "vcf"))]
  idx    <- x[["Qual"]] >= min_qual &
    x[["percentage_missing_samples"]] <= max_missing_samples &
    x[["MAF"]] >= maf_cutoff
  stopifnot("\nNo locus in VCF file has satisfied the specified QC metrics." =
              any(idx))
  if (sum(idx) < nrow(snpdata[["details"]])) {
    x    <- x[idx, ]
    snpdata[["details"]] <- x
    for (field in fields) {
      snpdata[[field]] <- snpdata[[field]][idx, ]
    }
  } else if (sum(idx) == nrow(snpdata[["details"]])) {
    message("All loci have satisfied the specified QC metrics.")
  }
  
  # filter the meta table
  idx <- snpdata[["meta"]][["percentage_missing_sites"]] <= max_missing_sites
  stopifnot("\nNo sample in VCF file has satisfied the specified QC metrics." =
              any(idx))
  if (sum(idx) < nrow(snpdata[["meta"]])) {
    samples_to_be_dropped <- snpdata[["meta"]][["sample"]][!idx]
    message("\nThe following samples will be removed:\n",
            paste(samples_to_be_dropped, collapse = ", "))
    snpdata[["meta"]] <- snpdata[["meta"]][idx, ]
  } else if (sum(idx) == nrow(snpdata[["meta"]])) {
    message("\nAll samples have satisfied the specified QC metrics.")
  }
  
  # remove bad quality samples from genotype matrices and recalculate both SNPs
  # and sample missingness and MAF

  # filter out samples with bad qualities from the genotype tables
  if (exists("samples_to_be_dropped")) {
    m <- match(samples_to_be_dropped, colnames(snpdata[["GT"]]))
    for (field in fields) {
      snpdata[[field]] <- snpdata[[field]][, -m]
    }

    # recalculate the percent of missing data for every SNP
    snpdata[["details"]][["percentage_missing_samples"]] <-
      rowSums(is.na(snpdata[["GT"]])) / ncol(snpdata[["GT"]])
    
    # recalculate the MAF
    snpdata <- compute_maf(
      snpdata,
      include_het = FALSE,
      mat_name    = "GT"
    )
  }
  
  # recalculate the percent of missing data per sample
  snpdata[["meta"]][["percentage_missing_sites"]] <-
    colSums(is.na(snpdata[["GT"]])) / nrow(snpdata[["GT"]])
  
  snpdata
}
