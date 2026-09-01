#' Remove a set of samples from a <SNPdata> object
#'
#' @param snpdata The input <SNPdata> object
#' @param samples A vector or a file with one column containing the list of
#'    samples to be dropped from the provided <SNPdata> object. When the value
#'    of this argument is a file, the column must be named as `sample`.
#'
#' @returns A <SNPdata> object where the specified samples have been removed
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
#'  # remove samples
#'  snpdata <- remove_samples(
#'    snpdata = snpdata,
#'    samples = c("PA0007-C", "PA0008-C")
#'  )
#' } 
remove_samples <- function(snpdata, samples) {
  checkmate::assert_class(snpdata, "SNPdata", null.ok = FALSE)
  if (!checkmate::check_vector(samples, any.missing = FALSE, null.ok = FALSE,
                               min.len = 1L)) {
    checkmate::assert_file_exists(samples, access = "r")
    target_samples <- data.table::fread(samples, header = TRUE)
    if (ncol(target_samples) != 1 && names(target_samples) != "sample") {
      cli::cli_abort(c(
        i = "Expected a file with one column named as {.val sample}",
        x = "Incorrect value for the {.strong samples} argument"
      ))
    }
    samples <- target_samples[["sample"]]
  }
  
  # check whether all provided samples are found in the 'meta' table
  if (!all(samples %in% snpdata[["meta"]][["sample"]])) {
    cli::cli_abort(c(
      i = "Some specified samples are not found in the input {.cls SNPdata} \\\
      object",
      x = "Invalid sample names were detected."
    ))
  }
  idx <- match(samples, snpdata[["meta"]][["sample"]])
  tmp_meta <- snpdata[["meta"]]
  tmp_meta <- tmp_meta[-(idx), ]
  snpdata[["meta"]] <- tmp_meta
  fields <- c("GT", "Phased", "Imputed")
  for (field in fields) {
    if (field %in% names(snpdata)) {
      idx <- match(samples, colnames(snpdata[[field]]))
      m <- seq_len(ncol(snpdata[[field]]))
      m <- m[-idx]
      tmp <- snpdata[[field]]
      tmp <- tmp[, m]
      snpdata[[field]] <- tmp
    }
  }
  
  # re-calculating the MAF from the 'GT' matrix
  if ("MAF" %in% names(snpdata[["details"]])) {
    snpdata <- compute_maf(
      snpdata = snpdata,
      genotype = "GT",
      include_het = FALSE
    )
  }
  
 # re-calculating the percent of missing genotypes for every SNPs from the 'GT'
 # matrix
 snpdata[["details"]][["percentage_missing_samples"]] <-
   rowSums(is.na(snpdata[["GT"]])) / ncol(snpdata[["GT"]])

  return(snpdata)
}