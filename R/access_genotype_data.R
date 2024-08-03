#' Extract the genotypes from a variant call file (VCF)
#'
#' @param file A character with the path to the input VCF file.
#'
#' @return A data frame with the genotype data together with the SNPs genomic
#'    coordinates and their associated calling quality
#' @export
#' @importFrom data.table %like%
#' @author Banky
#'
#' @examples
#' \dontrun{
#'   extract_genotype <- extract_allelic_depth(
#'     file = system.file("extdata", "Input_Data.vcf.gz", package = "mpbr")
#'   )
#' }
#'
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
  cli::cli_alert_success("\nThe sample genotypes were successfully extracted.")
  vcf
}


#' Extract the allelic depth from a variant call file (VCF)
#'
#' @param file A character with the path to the input VCF file.
#'
#' @return A data frame with the allelic depth supporting each of the two
#'    alleles at any given locus across all samples.
#' @export
#' @importFrom data.table %like%
#'
#' @examples
#' \dontrun{
#'   allelic_depth <- extract_allelic_depth(
#'     file = system.file("extdata", "Input_Data.vcf.gz", package = "mpbr")
#'   )
#' }
#'
extract_allelic_depth <- function(file) {
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
    , -c(1:9) # only keep the columns of the genotype field
  ]
  
  cols <- names(vcf)
  
  vcf[, (cols) := lapply(.SD, function(g) {
    unlist(strsplit(g, ":", fixed = TRUE))[[2L]]
  }), .SDcols = cols]
  cli::cli_alert_success("\nThe allelic depth were successfully extracted.")
  vcf
}