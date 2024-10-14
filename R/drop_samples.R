#' Drop a specific set of samples from the `SNPdata` object
#'
#' @param snpdata The input `SNPdata` object
#' @param samples_to_be_dropped A vector of samples to be dropped
#'
#' @return A `SNPdata` object where the specified samples have been removed
#' @examples
#' \dontrun{
#'   samples_to_be_dropped <- sample(
#'     snpdata[["meta"]][["sample"]],
#'     20,
#'     replace = FALSE
#'   )
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
  idx <- match(samples_to_be_dropped, snpdata[["meta"]][["sample"]])
  tmp_meta <- snpdata[["meta"]]
  tmp_meta <- tmp_meta[-(idx), ]
  snpdata[["meta"]] <- tmp_meta
  fields <- c("GT", "Phased", "Imputed")
  for (field in fields) {
    if (field %in% names(snpdata)) {
      idx <- match(samples_to_be_dropped, colnames(snpdata[[field]]))
      m <- seq_len(ncol(snpdata[[field]]))
      m <- m[-idx]
      tmp <- snpdata[[field]]
      tmp <- tmp[, m]
      snpdata[[field]] <- tmp
    }
  }
  
  message("Re-calculating the MAF from the 'GT' matrix...")
  snpdata <- compute_maf(
    snpdata = snpdata,
    include_het = FALSE,
    mat_name = "GT"
  )
  message("Re-calculating the percent of missing samples for every SNPs from",
          "the 'GT' matrix...")
  snpdata[["details"]][["percentage_missing_samples"]] <-
    rowSums(is.na(snpdata[["GT"]])) / ncol(snpdata[["GT"]])

  snpdata
}