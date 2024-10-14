#' Select data from specified chromosomes
#'
#' @param snpdata a `SNPdata` object
#' @param chrom a vector of chromosome names
#'
#' @return a `SNPdata` object with only the data from the specified chromosomes
#'
#' @examples
#' \dontrun{
#'   chrom7 <- select_chrom(snpdata, chrom = "Pf3D7_07_v3")
#'  }
#' @export
#'
select_chrom <- function(snpdata, chrom = "all") {
  checkmate::assert_class(snpdata, "SNPdata", null.ok = FALSE)
  checkmate::assert_vector(chrom, any.missing = FALSE, null.ok = FALSE,
                           min.len = 1L)
  # return the input object if no chromosome is specified
  if (chrom == "all") {
    message("No chromosome was specified. Returning the initial object.")
    return(snpdata)
  }

  # subset if a chromosome name is specified
  # the specified chromosomes must be valid names
  stopifnot("Specified chromosome name not found in the input SNPdata object"
            = chrom %in% unique(snpdata[["details"]][["Chrom"]]))
  res <- list()
  idx <- names(snpdata) %in% c("vcf", "index", "meta")
  fields <- names(snpdata)[!idx]
  for (chr in chrom) {
    chrom_snpdata <- snpdata
    idx <- which(chrom_snpdata[["details"]][["Chrom"]] == chr)
    for (field in fields) {
      res[[field]] <- rbind(res[[field]], chrom_snpdata[[field]][idx, ])
    }
  }

  # recalculate the percent of missing values per samples
  # this accounts for the drop in the number of SNPs
  meta <- snpdata[["meta"]]
  meta[["percentage_missing_sites"]] <-
    colSums(is.na(res[["GT"]])) / nrow(res[["GT"]])
  res[["meta"]] <- meta
  class(res) <- "SNPdata"
  res
}
