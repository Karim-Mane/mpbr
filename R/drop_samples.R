#' Drop a specific set of samples from the `SNPdata` object
#'
#' @param snpdata a `SNPdata` object
#' @param samples_to_be_dropped a vector of samples to be dropped
#'
#' @return a `SNPdata` object where the specified samples have been removed
#' @examples
#' \dontrun{
#'   snpdata <- drop_samples(snpdata, samples_to_be_dropped)
#'  }
#'
#' @export
drop_samples <- function(snpdata, samples_to_be_dropped) {
  checkmate::assert_class(snpdata, "SNPdata", null.ok = FALSE)
  checkmate::assert_vector(samples_to_be_dropped, any.missing = FALSE,
                           null.ok = FALSE, min.len = 1L)
  stopifnot("Some samples in the provided vector are not found in the current
         data" = all(samples_to_be_dropped %in% snpdata[["meta"]][["sample"]]))
  idx               <- match(samples_to_be_dropped,
                             snpdata[["meta"]][["sample"]])
  tmp_meta          <- snpdata[["meta"]]
  tmp_meta          <- tmp_meta[-(idx), ]
  snpdata[["meta"]] <- tmp_meta
  fields            <- c("GT", "Phased", "Imputed")
  for (field in fields) {
    if (field %in% names(snpdata)) {
      idx              <- match(samples_to_be_dropped,
                                colnames(snpdata[[field]]))
      m                <- seq_len(ncol(snpdata[[field]]))
      m                <- m[-idx]
      tmp              <- snpdata[[field]]
      tmp              <- tmp[, m]
      snpdata[[field]] <- tmp
    }
  }
  tmp_file <- file.path(dirname(snpdata[["vcf"]]), "tmp.txt")
  data.table::fwrite(snpdata[["meta"]][["sample"]], tmp_file, col.names = FALSE,
                     row.names = FALSE, quote = FALSE, sep = "\t")
  index    <- ifelse("index" %in% names(snpdata),
                     snpdata[["index"]] + 1L,
                     0L)
  snpdata[["vcf"]]  <- remove_samples_from_vcf(snpdata[["vcf"]], "tmp.txt",
                                               path = dirname(snpdata[["vcf"]]),
                                               index = index)
  snpdata[["index"]] <- index
  snpdata
}

#' Remove samples from VCF file
#'
#' @param vcf the input VCF file
#' @param samples_to_be_retained a vector of samples to be retained
#' @param path the path folder that contains the input data. This will be also
#'    the folder where the output files will be stored.
#' @param index the index of the input VCF file.
#'
#' @return the path to filtered VCF file
#' @keywords internal
#'
remove_samples_from_vcf <- function(vcf, samples_to_be_retained,
                                    path, index = 1L) {
  checkmate::assert_file_exists(vcf)
  checkmate::assert_vector(samples_to_be_retained, any.missing = FALSE,
                           null.ok = FALSE, min.len = 1L)
  checkmate::assert_directory_exists(path)
  checkmate::assert_numeric(index, lower = 0L, any.missing = FALSE,
                            null.ok = FALSE, len = 1L)
  target_samples <- file.path(path, samples_to_be_retained)
  post_qc        <- file.path(path, paste0("Post_QC_", index, ".vcf.gz"))
  system(sprintf("bcftools view -S %s %s -o %s -O z",
                 target_samples, vcf, post_qc))
  system(sprintf("tabix %s", post_qc))
  as.character(post_qc)
}