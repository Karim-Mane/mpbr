#' Extract either the genotypes (GT) or the allelic depth (AD) from a variant
#' call file (VCF).
#'
#' @param vcf_data A <data table> object obtained after reading the input VCF
#'    file.
#' @param which A character that denotes the data of interest. The possible
#'    values are:
#'   \enumerate{
#'      \item GT: to extract the genotype data (FORMAT/GT)
#'      \item AD: to extract the allelic depth (FORMAT/AD)
#'   }
#' @param nthreads A numeric that represents the number of threads to use during
#'    data import and processing
#'
#' @return A data frame with the genotype data of interest, together with the
#'    SNPs genomic coordinates and their associated alleles (if which = "GT").
#' @export
#' @importFrom data.table %like%
#' @author Banky
#' @examples
#' \dontrun{
#'   # read in the VCF file
#'   vcf_data <- data.table::fread(
#'     file = system.file("extdata", "test_data.vcf.gz", package = "mpbr"),
#'     skip = "#CHROM",
#'     header = TRUE,
#'     sep = "\t",
#'     nThread = nthreads
#'   )
#'   # extract genotype matrix
#'   genotypes <- extract_genotype(
#'     vcf_data = vcf_data,
#'     which = "GT",
#'     nthreads = 5
#'   )
#'   # extract allelic depth matrix
#'   genotypes <- extract_genotype(
#'     vcf_data = vcf_data,
#'     which = "AD",
#'     nthreads = 5
#'   )
#' }
extract_genotype <- function(vcf_data, which = "GT", nthreads) {
  checkmate::assert_data_table(vcf_data, null.ok = FALSE)
  checkmate::assert_character(which, n.chars = 2L, len = 1L, null.ok = FALSE,
                              any.missing = FALSE)
  which <- match.arg(which, choices = c("GT", "AD"))

  # removing the unnecessary columns. When 'which = "GT"', we don't want to keep
  # the ID, QUAL, FILTER, INFO, and FORMAT columns. For "AD", only the CHROM,
  # POS, and ADs across every sample will be returned. This avoids extracting
  # redundant information which would have already been collected in the GT
  # data.
  idx <- c(3, 6, 7, 8, 9)
  idx_genotype_data <- 5
  genotype <- "genotypes"
  if (which == "AD") {
    idx <- 3:9
    idx_genotype_data <- 3
  }
  vcf_data <- subset(vcf_data, select = -idx)
  
  # element of the FORMAT/GT field are separated by ':'. The genotype (GT) field
  # is the first element, but the allelic depth (AD) is the second element
  get_gt_field <- function(x) {
    x <- as.matrix(x)
    sub(":.*", "", x, perl = TRUE) # keep text before first ":"
  }
  genotypes <- get_gt_field(vcf_data[, idx_genotype_data:length(vcf_data)])
  genotypes <- cbind(vcf_data[, 1:(idx_genotype_data - 1)], genotypes)
  return(as.data.frame(genotypes))
}


# Extract SNPs functional effect annotations from INFO/ANN field of a VCF file.
# These annotations will be added to the 'details' table of the <SNPdata>
# object.
#'
#' @inheritParams extract_genotype vcf_data nthreads
#'
#' @returns A data frame with eight columns including each variant annotation
#'    effect
#' @export
#'
#' @examples
#' \dontrun{
#'   # read in the VCF file
#'   vcf_data <- data.table::fread(
#'     file = system.file("extdata", "test_data.vcf.gz", package = "mpbr"),
#'     skip = "#CHROM",
#'     header = TRUE,
#'     sep = "\t",
#'     nThread = nthreads
#'   )
#'   # rename #CHROM to CHROM
#'   data.table::setnames(vcf_data, "#CHROM", "CHROM")
#'   # extract annotations
#'   snp_effect_annot <- extract_annotation(
#'     vcf_data = vcf_data,
#'     nthreads = 5
#'   )
#' }
extract_annotation <- function(vcf_data, nthreads) {
  # extract ANN from INFO column
  # INFO looks like: ...;ANN=allele|effect|impact|...,allele|effect|...;...
  vcf_data[, ann_raw := sub(".*ANN=([^;]+).*", "\\1", INFO)]
  
  # parse all annotations per site
  ann_dt <- vcf_data[, {
    # split multiple annotations (comma-separated)
    entries <- unlist(strsplit(ann_raw, ",", fixed = TRUE))
    # split each annotation by pipe and extract fields 2,3,10,11
    fields <- lapply(entries, function(e) {
      f <- strsplit(e, "|", fixed = TRUE)[[1]]
      list(
        effect = if (length(f) >= 2)  f[2]  else NA_character_,
        impact = if (length(f) >= 3)  f[3]  else NA_character_,
        hgvs_c = if (length(f) >= 10) f[10] else NA_character_,
        hgvs_p = if (length(f) >= 11) f[11] else NA_character_
      )
    })
    # collapse unique values into comma-separated strings
    list(
      effect = paste(
        unique(na.omit(sapply(fields, `[[`, "effect"))),
        collapse = ","
      ),
      impact = paste(
        unique(na.omit(sapply(fields, `[[`, "impact"))),
        collapse = ","
      ),
      hgvs_c = paste(
        unique(na.omit(sapply(fields, `[[`, "hgvs_c"))),
        collapse = ","
      ),
      hgvs_p = paste(
        unique(na.omit(sapply(fields, `[[`, "hgvs_p"))),
        collapse = ","
      )
    )
  }, by = .(CHROM, POS, REF, ALT)]

  # clean empty strings and dots
  ann_dt[ann_dt == "" | ann_dt == "."] <- NA

  return(ann_dt)
}