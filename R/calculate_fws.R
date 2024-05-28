#' Calculate the complexity of the infection (COI) in every sample
#'
#' The COI is estimated here through the within host genetic diversity (Fws).
#' using the {moimix} R package.
#'
#' @param snpdata a `SNPdata` object
#' @param threshold a `numeric` value between 0 and 1. This will be used to
#'    categorize the infections into monoclonal (Fws > threshold) or polyclonal
#'    (Fws <= threshold)
#'
#' @return a `SNPdata` object with 2 additional columns in the meta table
#' \enumerate{
#'   \item Fws: within host genetic diversity value
#'   \item COI: the complexity of infection: 1 for Fws>0.95, 2 for Fws<=0.95
#' }
#'
#' @examples
#' \dontrun{
#'   snpdata <- calculate_fws(snpdata, threshold = 0.95)
#'  }
#'
#' @export
#'
calculate_fws <- function(snpdata, threshold = 0.95) {
  checkmate::assert_class(snpdata, "SNPdata", null.ok = FALSE)
  checkmate::assert_numeric(threshold, lower = 0L, upper = 1L, finite = TRUE,
                            len = 1L, null.ok = FALSE)
  vcf        <- snpdata[["vcf"]]
  gds_file   <- file.path(dirname(vcf), "data.gds")
  SeqArray::seqVCF2GDS(vcf, gds_file)
  my_vcf    <- SeqArray::seqOpen(gds_file)
  
  # filter out the SNPs from the apicoplast and the mitochondrial chromosomes
  sample_id <- SeqArray::seqGetData(my_vcf, "sample.id")
  coords    <- moimix::getCoordinates(my_vcf)
  SeqArray::seqSetFilter(my_vcf,
                         variant.id = coords[["variant.id"]][coords[["chromosome"]] != "Pf3D7_API_v3"]) # nolint: line_length_linter
  SeqArray::seqSetFilter(my_vcf,
                         variant.id = coords[["variant.id"]][coords[["chromosome"]] != "Pf_M76611"]) # nolint: line_length_linter
  
  # calculate the Fws and the COI
  fws_overall        <- moimix::getFws(my_vcf)
  fws_overall        <- data.table::data.table(cbind(as.character(sample_id),
                                                     as.numeric(fws_overall)))
  names(fws_overall) <- c("sample", "Fws")
  meta <- data.frame(snpdata[["meta"]]) %>%
    dplyr::left_join(fws_overall, by = "sample")
  meta[["COI"]]      <- 1L
  meta[which(meta[["Fws"]] <= threshold), ][["COI"]] <- 2L
  snpdata[["meta"]]  <- meta
  snpdata
}