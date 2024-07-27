#' Extract the genotypes from a variant call file (VCF)
#'
#' @param file A character with the path to the input VCF file.
#'
#' @return A data frame with the genotype data together with the SNPs genomic
#'    coordinates and their associated calling quality
#' @keywords internal
#' @importFrom data.table %like%
#' @author Banky
extract_genotype <- function(file) {
  vcf <- data.table::fread(
    cmd     = sprintf("pigz -dc < %s", file),
    nThread = (parallel::detectCores() - 2),
    sep     = NULL, # line is column
    header  = FALSE # 1st line is data
  )
  
  # DT name for unnamed col
  V1 <- NULL # nolint: object_name_linter.
  row_id <- vcf[V1 %like% "^#CHROM", which = TRUE] + 1 # + 1 skip header
  
  vcf <- vcf[row_id:.N, data.table::tstrsplit(V1, "\t", fixed = TRUE)][
    , -c(3, 7, 8, 9) # remove unwanted columns ID, FILTER, INFO, FORMAT
  ]
  
  cols <- names(vcf)[6:length(vcf)] # first 6 don't need processing
  
  vcf[, (cols) := lapply(.SD, function(g) {
    substring(g, 1, regexpr(":", g, fixed = TRUE) - 1) # 1st position only
  }), .SDcols = cols]
  vcf
}