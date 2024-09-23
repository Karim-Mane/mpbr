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
  ## check your operating system type
  os_type <- Sys.info()["sysname"]
  
  if (os_type == "Windows") {
    vcf <- data.table::fread(
      cmd     = sprintf("gzip -dc < %s", file),
      nThread = (parallel::detectCores() - 2),
      sep     = NULL, # line is column
      header  = FALSE # 1st line is data
    )
  } else {
    vcf <- data.table::fread(
      cmd     = sprintf("pigz -dc < %s", file),
      nThread = (parallel::detectCores() - 2),
      sep     = NULL, # line is column
      header  = FALSE # 1st line is data
    )
  }
  
  # find the row number where the header #CHROM Starts.
  row_id <- which(grepl("^#CHROM", vcf$V1)) + 1
  
  # splot the single column V1 into multiple columns then later delet unwanted columns
  vcf <- vcf[row_id:.N, data.table::tstrsplit(V1, "\t", fixed = TRUE)][
    , -c(3, 7, 8, 9) # remove unwanted columns ID, FILTER, INFO, FORMAT
  ]
  
  # select the columns that contain the genotype info
  cols <- names(vcf)[6:length(vcf)] # first 5 don't need processing
  
  # extract the genotype info from th slected columns
  vcf[, (cols) := lapply(.SD, function(g) {
    substring(g, 1, regexpr(":", g, fixed = TRUE) - 1) # 1st position only
  }), .SDcols = cols]
  vcf
}