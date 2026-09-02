#' Filter a <SNPdata> object
#'
#' SNPs and samples are filtered based on the user-specified thresholds. Both
#' the `meta`, `details` and the different genotype tables will be affected
#' depending on the filtration conditions.
#'
#' @param snpdata An object of class <SNPdata>
#' @param max_missing_sites A numeric representing the maximum fraction of
#'    missing sites above which a sample should be dropped. Default is 0.2.
#' @param max_missing_samples A numeric representing the maximum fraction of
#'    missing samples above which a loci should be discarded. Default is 0.2.
#' @param maf_cutoff A numeric representing the minor allele frequency cut-off.
#'    Loci with a `MAF < maf_cutoff` will be discarded.
#'
#' @return a filtered <SNPdata> object
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
#'  # filter the SNPdata object
#'  snpdata <- filter(
#'    snpdata = snpdata,
#'    max_missing_sites = 0.2,
#'    max_missing_samples = 0.2,
#'    maf_cutoff = 0.01
#'  )
#' }
filter <- function(snpdata,
                   max_missing_sites = 0.2,
                   max_missing_samples = 0.2,
                   maf_cutoff = 0.01) {
  checkmate::assert_class(snpdata, "SNPdata", null.ok = FALSE)
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
  fields <- names(snpdata)[!(names(snpdata) %in% c("meta", "details", "vcf"))]
  # calculate the MAF if it was not
  if (!"MAF" %in% names(snpdata[["details"]])) {
    snpdata <- compute_maf(
      snpdata = snpdata,
      include_het = FALSE,
      genotype = "GT"
    )
  }
  x <- snpdata[["details"]]
  idx <- x[["percentage_missing_samples"]] <= max_missing_samples &
    x[["MAF"]] >= maf_cutoff
  if (sum(idx) == 0) {
    cli::cli_abort(c(
      i = "No SNP in the input <SNPdata> object has satisfied the specified \\\
      filtering thresholds.",
      x = "SNPs filtering thresholds might be too high. Consider using less \\\
      strict thresholds."
    ))
  }
  if (sum(idx) < nrow(snpdata[["details"]])) {
    x <- x[idx, ]
    snpdata[["details"]] <- x
    for (field in fields) {
      snpdata[[field]] <- snpdata[[field]][idx, ]
    }
  } else if (sum(idx) == nrow(snpdata[["details"]])) {
    cli::cli_alert_success(
    "All SNPs have satisfied the specified filtering thresholds")
  }
  
  # filter the meta table
  idx <- snpdata[["meta"]][["percentage_missing_sites"]] <= max_missing_sites
  if (sum(idx) == 0) {
    cli::cli_abort(c(
      i = "No sample in input <SNPdata> object has satisfied the specified \\\
      filtering thresholds.",
      x = "samples filtering thresholds might be too high. Consider using \\\
      less strict thresholds."
    ))
  }
  samples_to_be_dropped <- character()
  if (sum(idx) < nrow(snpdata[["meta"]])) {
    samples_to_be_dropped <- snpdata[["meta"]][["sample"]][!idx]
    cli::cli_inform(
      "The following samples will be removed: \\\
      {.val {toString(samples_to_be_dropped)}}",
    )
    snpdata[["meta"]] <- snpdata[["meta"]][idx, ]
  } else if (sum(idx) == nrow(snpdata[["meta"]])) {
    cli::cli_alert_success(
      "All samples have satisfied the specified filtering thresholds."
    )
  }
  
  # remove bad quality samples from genotype matrices and recalculate both SNPs
  # and sample missingness and MAF

  # filter out samples with bad qualities from the genotype tables
  if (length(samples_to_be_dropped) > 0) {
    m <- match(samples_to_be_dropped, colnames(snpdata[["GT"]]))
    for (field in fields) {
      snpdata[[field]] <- snpdata[[field]][, -m]
    }

    # recalculate the percent of missing genotypes for every SNP
    snpdata[["details"]][["percentage_missing_samples"]] <-
      rowSums(is.na(snpdata[["GT"]])) / ncol(snpdata[["GT"]])
    
    # recalculate the MAF
    snpdata <- compute_maf(
      snpdata,
      include_het = FALSE,
      genotype    = "GT"
    )
  }
  
  # recalculate the percent of missing data per sample
  snpdata[["meta"]][["percentage_missing_sites"]] <-
    colSums(is.na(snpdata[["GT"]])) / nrow(snpdata[["GT"]])
  
  return(snpdata)
}
