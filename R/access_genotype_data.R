#' Extract either the genotypes (GT) or the allelic depth (AD) from a variant
#' call file (VCF).
#'
#' @param file A character with the path to the input VCF file.
#' @param which A character that denotes the data of interest. The possible
#'    values are:
#'   \enumerate{
#'      \item GT: to extract the genotype data
#'      \item AD: to extract the allelic depth
#'   }
#'
#' @return A data frame with the genotype data of interest, together with the
#'    SNPs genomic coordinates and their associated alleles and calling quality
#'    (if which = "GT").
#' @export
#' @importFrom data.table %like%
#' @author Banky
#'
extract_genotype <- function(file, which = "GT") {
  checkmate::assert_file_exists(file)
  checkmate::assert_character(which, n.chars = 2L, len = 1L, null.ok = FALSE,
                              any.missing = FALSE)
  which <- match.arg(which, choices = c("GT", "AD"))

  # import the VCF file using the {data.table} package and pigz.
  # users who have enabled openmp. will benefit from {data.table} multithreaded
  # data import 
  vcf <- data.table::fread(
    cmd = sprintf("pigz -dc < %s", file),
    nThread = (parallel::detectCores() - 2),
    sep = NULL, # line is column
    header = FALSE # 1st line is data
  )

  # there are 2 parts in the vcf file: the header (rows with '##' or '#') and
  # the genotype (rows after the headers). The below extract the row number
  # (row_id) from where the genotype section starts. The data of interest is
  # located between that row and the last row in the file.
  V1 <- NULL # nolint: object_name_linter.
  row_id <- vcf[V1 %like% "^#CHROM", which = TRUE] + 1 # + 1 skip header

  # removing the unnecessary columns. When 'which = "GT"', we don't want to keep
  # the ID, FILTER, INFO, and FORMAT columns. For "AD", only the CHROM, POS, and
  # ADs across every sample will be returned. This avoids extracting redundant
  # information which would have already been collected in the GT data.
  idx <- c(3, 7, 8, 9)
  idx_genotype_data <- 6
  idx_extraction <- 1
  genotype <- "genotypes"
  if (which == "AD") {
    idx <- 3:9
    idx_genotype_data <- 3
    idx_extraction <- 2
    genotype <- "allelic depths"
  }
  vcf <- vcf[row_id:.N, data.table::tstrsplit(V1, "\t", fixed = TRUE)]
  vcf <- subset(vcf, select = -idx)

  # The GT data is located from the 6th column of the data, while the AD is
  # found from column 3. See the section above for the reasoning.
  genotypes <- names(vcf)[idx_genotype_data:length(vcf)]
  
  # the genotype field is split based on ':' as a pattern. The genotypes data
  # is the first element, but the allelic depth is the second element
  get_gt_field <- function(x) {
    x <- as.matrix(x)
    sub(":.*", "", x, perl = TRUE) # keep text before first ":"
  }
  genotypes <- get_gt_field(vcf[, idx_genotype_data:length(vcf)])
  genotypes <- cbind(vcf[, 1:(idx_genotype_data - 1)], genotypes)
  return(as.data.frame(genotypes))
}