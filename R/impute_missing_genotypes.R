#' Impute missing genotypes
#'
#' Missing genotype imputation is performed based on the MAF at any given locus.
#' Imputation will be done `nsim` times and imputed data with highest
#' correlation between MAF from raw data and MAF from imputed data will be
#' retained
#'
#' @param snpdata a `SNPdata` object
#' @param genotype the genotype table from which the missing data will
#'    be imputed. This can be either the raw genotype matrix (`GT`) or the
#'    phased genotype matrix (`Phased`)
#' @param nsim an integer that represents the number of simulations
#'
#' @return a `SNPdata` object with an additional table named as:
#'    "Phased_Imputed" if the phased data was used for imputation or "Imputed"
#'    if the imputation was done on the raw genotypes.
#'
#' @details When both alleles are not supported by any read or the total number
#'    of reads supporting both alleles at a given site is < 5, the genotype will
#'    be phased based on a Bernoulli distribution using the MAF as a parameter.
#'    Similarly, when the total number of reads is > 5 and the number of reads
#'    supporting one of the allele is not 2 times the number of the other
#'    allele, the genotype is phased using a Bernoulli distribution.
#'
#' @examples
#' \dontrun{
#'   snpdata <- impute_missing_genotypes(snpdata)
#'  }
#'
#' @export
#'
impute_missing_genotypes <- function(snpdata, genotype = "Phased",
                                     nsim = 100L) {
  checkmate::assert_class(snpdata, "SNPdata", null.ok = FALSE)
  checkmate::assert_character(genotype, any.missing = FALSE, len = 1L,
                              null.ok = FALSE)
  checkmate::assert_numeric(nsim, lower = 1L, any.missing = FALSE,
                            null.ok = FALSE, len = 1L)
  message("The missing genotypes will be imputed from ", genotype, " table.\n")
  field <- genotype
  path  <- file.path(dirname(snpdata[["vcf"]]), "imputing")
  system(sprintf("mkdir -p %s", path))
  correlations <- numeric(length = nsim)
  pb           <- utils::txtProgressBar(min = 0L, max = nsim, initial = 0L,
                                        style = 3L, char = "*")
  for (i in seq_len(nsim)) {
    tmp_snpdata <- snpdata
    mat         <- apply(tmp_snpdata[[field]], 1L, impute)
    tmp_snpdata[["Imputed"]] <- t(mat)
    saveRDS(t(mat), file.path(path, paste0("sim", i, ".RDS")))
    res_snpdata <- compute_maf(tmp_snpdata, include_het = FALSE,
                               mat_name = "Imputed")
    correlations[i] <- stats::cor(res_snpdata[["details"]][["MAF_Imputed"]],
                                  res_snpdata[["details"]][["MAF"]])
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
  idx <- which(correlations == max(correlations, na.rm = TRUE))
  snpdata[["Imputed"]] <- readRDS(file.path(path,
                                            paste0("sim", idx[[1L]], ".RDS")))
  system(sprintf("rm -rf %s", path))
  snpdata
}

#' Impute missing genotypes for 1 SNP across all samples
#'
#' @param genotype a vector of intergers
#'
#' @return the input vector where all missing alleles have been imputed.
#' @keywords internal
#'
impute <- function(genotype) {
  checkmate::assert_vector(genotype, min.len = 1L, null.ok = FALSE)
  idx <- as.numeric(which(is.na(genotype)))
  for (j in idx) {
    ref         <- length(which(genotype == 0L))
    alt         <- length(which(genotype == 1L))
    maf         <- ifelse(ref < alt, ref / (ref + alt), alt / (ref + alt))
    genotype[j] <- statip::rbern(1L, maf)
  }
  genotype
}