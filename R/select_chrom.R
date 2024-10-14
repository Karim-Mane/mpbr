#' Select data from specified chromosomes
#'
#' @param snpdata a `SNPdata` object
#' @param chrom a vector of chromosome names
#'
#' @return a `SNPdata` object with only the data from the specified chromosomes
#'
#' @examples
#' \dontrun{
#'   chrom_snpdata <- select_chrom(snpdata, chrom = "Pf3D7_07_v3")
#'  }
#' @export
#'
select_chrom <- function(snpdata, chrom = "all") {
  checkmate::assert_class(snpdata, "SNPdata", null.ok = FALSE)
  checkmate::assert_vector(chrom, any.missing = FALSE, null.ok = FALSE,
                           min.len = 1L)
  system(sprintf("tabix -f %s", snpdata[["vcf"]]))
  ## make sure to output the metadata
  percentage_missing_sites <- Fws <- COI <- NULL
  meta   <- snpdata[["meta"]] %>%
    dplyr::select(-c(percentage_missing_sites, Fws, COI)) # nolint: object_name_linter

  # return the input object if no chromosome is specified
  if (chrom == "all") {
    return(snpdata)
  }

  # subset if a chromosome name is specified
  stopifnot("Specified chromosome name not found in the input SNPdata object"
            = chrom %in% unique(snpdata[["details"]][["Chrom"]]))
  res <- list()
  m <- which(names(snpdata) %in% c("vcf", "index", "meta"))
  fields <- names(snpdata)[-m]
  for (chr in chrom) {
    chrom_snpdata <- snpdata
    idx <- which(chrom_snpdata[["details"]][["Chrom"]] == chr)
    for (field in fields) {
      res[[field]] <- rbind(res[[field]], chrom_snpdata[[field]][idx, ])
    }
  }
  chrom_vcf <- file.path(dirname(chrom_snpdata[["vcf"]]), "target_chrom_vcf.gz")
  if (file.exists(chrom_vcf)) {
    unlink(chrom_vcf)
  }
  if (length(chrom) > 1L) {
    tmp_xme <- file.path(dirname(chrom_snpdata[["vcf"]]), "target_chrom.txt")
    data.table::fwrite(chrom, tmp_xme, col.names = FALSE,
                       row.names = FALSE, quote = FALSE, sep = "\t")
    system(sprintf("bcftools view -R %s %s -o %s -O z",
                   tmp_xme, chrom_snpdata[["vcf"]], chrom_vcf))
  } else {
    system(sprintf("bcftools view -r\"%s\" %s -o %s -O z",
                   chrom, chrom_snpdata[["vcf"]], chrom_vcf))
  }
  res[["vcf"]] <- chrom_vcf
  res[["meta"]] <- meta
  class(res) <- "SNPdata"
  res
}
