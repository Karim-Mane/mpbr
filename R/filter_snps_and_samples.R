#' Filter loci and samples (requires **bcftools** and **tabix** to be installed)
#'
#' This function filters the SNPs and samples based on the specified conditions
#'
#' @param snpdata a `SNPdata` object
#' @param min_qual the minimum call quality score below which a loci will be
#'    discarded. default = 10
#' @param max_missing_sites the maximum fraction of missing sites above which
#'    a sample should be discarded. default = 0.2
#' @param max_missing_samples the maximum fraction of missing samples above
#'    which a loci should be discarded. default = 0.2
#'  @param iterative_fileration enables recursive filtration of missingness
#'   for loci and samples. default = TRUE 
#' @param maf_cutoff the MAF cut-off. loci with a MAF < maf_cutoff will be
#'    discarded
#'
#' @return a filtered SNPdata object
#' @examples
#' \dontrun{
#'  snpdata <- filter_snps_samples(
#'   snpdata,
#'   min_qual            = 10,
#'   max_missing_sites   = 0.2,
#'   max_missing_samples = 0.2,
#'   iterative_fileration = TRUE,
#'   maf_cutoff          = 0.01
#'  )
#'  }
#' @export
#'
filter_snps_samples <- function(snpdata,
                                min_qual            = 10L,
                                max_missing_sites   = 0.2,
                                max_missing_samples = 0.2,
                                iterative_fileration = TRUE,
                                maf_cutoff          = 0.01) {
  checkmate::assert_class(snpdata, "SNPdata", null.ok = FALSE)
  checkmate::assert_numeric(min_qual, lower = 10L, any.missing = FALSE,
                            null.ok = FALSE)
  checkmate::assert_numeric(max_missing_sites, lower = 0L, upper = 1L,
                            finite = TRUE, any.missing = FALSE, null.ok = FALSE,
                            len = 1L)
  checkmate::assert_numeric(max_missing_samples, lower = 0L, upper = 1L,
                            finite = TRUE, any.missing = FALSE, null.ok = FALSE,
                            len = 1L)
  checkmate::assert_logical(iterative_fileration, any.missing = FALSE,
                            null.ok = FALSE)
  checkmate::assert_numeric(maf_cutoff, lower = 0L, upper = 1L,
                            finite = TRUE, any.missing = FALSE, null.ok = FALSE,
                            len = 1L)
  
  if(iterative_fileration == TRUE){
    x <- snpdata[["details"]]
    y <- snpdata[["meta"]]
    snps <- snpdata[["GT"]]
    if (max_missing_samples <= 80)
    {
      snps <- snps[which(x[["percentage_missing_samples"]] <= 0.8), which(y[["percentage_missing_sites"]] <= 0.8)]
      x= x[x[["percentage_missing_samples"]] <= 0.8,]
      y= y[y[["percentage_missing_sites"]] <= 0.8,]
      #compute missingness
      y[["percentage_missing_sites"]]      <- colSums(is.na(snps)) / nrow(snps)
      x[["percentage_missing_samples"]] <- rowSums(is.na(snps)) / ncol(snps)
    }
    
    if (max_missing_samples <= 50)
    {
      #filter 50% missingness for sample and loci
      snps <- snp[which(x[["percentage_missing_samples"]] <= 0.5), which(y[["percentage_missing_sites"]] <= 0.5)]
      x= x[x[["percentage_missing_samples"]] <= 0.5,]
      y= y[y[["percentage_missing_sites"]] <= 0.5,]
      #compute missingness
      y[["percentage_missing_sites"]]      <- colSums(is.na(snps)) / nrow(snps)
      x[["percentage_missing_samples"]] <- rowSums(is.na(snps)) / ncol(snps)
      snpdata[["details"]] <- x
      snpdata[["meta"]] <- x
      snpdata[["GT"]] <- snps
      
    }
  }
  else{
    x      <- snpdata[["details"]]
    
  }
  
  fields <- c("GT", "Phased", "Phased_Imputed")
  idx    <- which(x[["Qual"]] >= min_qual &
                    x[["percentage_missing_samples"]] <= max_missing_samples &
                    x[["MAF"]] >= maf_cutoff)
  stopifnot("\nNo locus in VCF file has satisfied the specified QC metrics." =
              length(idx) != 0L)
  if (all(length(idx) > 0L && length(idx) < nrow(snpdata[["details"]]))) {
    x    <- x[idx, ]
    snpdata[["details"]] <- x
    for (field in fields) {
      if (field %in% names(snpdata)) {
        snpdata[[field]] <- snpdata[[field]][idx, ]
      }
    }
    Chrom <- Pos <- NULL # nolint: object_name_linter
    f2c          <- x %>% dplyr::select(Chrom, Pos)
    output_dir   <- dirname(snpdata[["vcf"]])
    data.table::fwrite(f2c, file.path(output_dir, "loci_to_be_retained.txt"),
                       col.names = FALSE,
                       row.names = FALSE,
                       quote     = FALSE,
                       sep       = "\t",
                       nThread   = 4L)
    snpdata[["vcf"]] <- remove_snps_from_vcf(snpdata[["vcf"]],
                                             "loci_to_be_retained.txt",
                                             output_dir,
                                             index = snpdata[["index"]])
  } else if (length(idx) == nrow(snpdata[["details"]])) {
    message("All loci have satisfied the specified QC metrics.")
  }
  
  idx <- which(snpdata[["meta"]][["percentage_missing_sites"]] <=
                 max_missing_sites)
  stopifnot("\nNo sample in VCF file has satisfied the specified QC metrics." =
              length(idx) != 0L)
  if (all(length(idx) > 0L && length(idx) < nrow(snpdata[["meta"]]))) {
    message("\nthe following samples will be removed:\n",
            paste(snpdata[["meta"]][["sample"]][-idx], collapse = "\n"))
    snpdata[["meta"]] <- snpdata[["meta"]][idx, ]
    data.table::fwrite(snpdata[["meta"]][["sample"]][-idx],
                       file.path(output_dir, "samples_to_be_dropped.txt"),
                       col.names = FALSE, row.names = FALSE,
                       quote = FALSE, sep = "\t", nThread = 4L)
    snpdata[["vcf"]] <- remove_samples_from_vcf(snpdata[["vcf"]],
                                                "samples_to_be_dropped.txt",
                                                output_dir,
                                                index = snpdata[["index"]])
  } else if (length(idx) == nrow(snpdata[["meta"]])) {
    message("\nAll samples have satisfied the specified QC metrics.")
  }
  
  snpdata[["index"]] <- snpdata[["index"]] + 1L
  snpdata
}
