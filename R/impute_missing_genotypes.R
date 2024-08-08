#' Impute missing genotypes
#'
#' Missing genotype imputation is performed based on the MAF at any given locus.
#' Imputation will be done `nsim` times and imputed data with highest
#' correlation between MAF from raw data and MAF from imputed data will be
#' retained
#'
#' @param snpdata The input `SNPdata` object.
#' @param genotype The name of the genotype table from which the missing data
#'    will be imputed. This can be either the raw genotype matrix (`GT`) or the
#'    phased genotype matrix (`Phased`).
#' @param nsim An integer that represents the number of simulations.
#' @param maf_col_name A character with the name of the column with the MAF in
#'    the 'details' table. Default is 'MAF'.
#'
#' @return A `SNPdata` object with an additional table named as:
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
#' @export
#'
impute <- function(snpdata, genotype = "GT", nsim = 100L,
                   maf_col_name = "MAF") {
  checkmate::assert_class(snpdata, "SNPdata", null.ok = FALSE)
  checkmate::assert_character(genotype, any.missing = FALSE, len = 1L,
                              null.ok = FALSE)
  checkmate::assert_numeric(nsim, lower = 1L, any.missing = FALSE,
                            null.ok = FALSE, len = 1L)
  # Do not proceed if the specified genotype matrix does not exist
  if (!(genotype %in% names(snpdata))) {
    current_gt_matrices <- names(snpdata)[!(names(snpdata) %in% # nolint: object_usage_linter
                                              c("meta", "details", "vcf"))]
    cli::cli_abort(
      c("x" = "The specified genotype matrix {.code genotype} does not exist",
        "i" = "Current genotype matrices are: {.code {current_gt_matrices}}")
    )
  }

  # Do not proceed if the chosen genotype matrix does not contain any missing
  # data.
  if (!anyNA(snpdata[[genotype]])) {
    cli::cli_abort(
      c("x" = "No missing data found in specified genotype matrix: {genotype}",
        "i" = "")
    )
  }

  # Do not proceed if the details table of the snpdata object does not have a
  # MAF column
  if (!(maf_col_name %in% names(snpdata[["details"]]))) {
    cli::cli_abort(
      c("x" = "{maf_col_name} column not found in {.code snpdata$details}.",
        "i" = "Use the {.code compute_maf()} function to calculate the MAF or
        specify the name of the column with the MAF using the
        {.code maf_col_name} argument.")
    )
  }

  # create a temporary directory to store temporary files
  temp_dir <- tempdir()
  path     <- file.path(temp_dir, "imputing")
  dir.create(path)

  # The correlations vector will store the correlation coefficient between the
  # initial MAF and the MAF after every simulation.
  # The imputed data from the simulation where the MAF is the closest to the
  # initial MAF will be chosen. If several of them have the same MAF, the first
  # will be retained.
  correlations <- numeric(length = nsim)
  for (i in cli::cli_progress_along(seq_len(nsim),
                                    "Imputing missing genotypes")) {
    tmp_snpdata <- snpdata
    mat         <- t(apply(tmp_snpdata[[genotype]], 1L, impute_one_allele))
    tmp_snpdata[["Imputed"]] <- mat
    saveRDS(mat, file.path(path, paste0("sim", i, ".RDS")))
    res_snpdata <- compute_maf(tmp_snpdata,
                               include_het = FALSE,
                               genotype    = "Imputed",
                               name        = "MAF_Imputed")
    correlations[i] <- stats::cor(res_snpdata[["details"]][["MAF_Imputed"]],
                                  res_snpdata[["details"]][["MAF"]])
  }
  cli::cli_alert_success("\nMissing genotypes were successfully imputed.")
  idx <- which(correlations == max(correlations, na.rm = TRUE))
  snpdata[["Imputed"]] <- readRDS(
    file.path(path, paste0("sim", idx[[1L]], ".RDS"))
  )

  # delete the temporary directory and return the output snpdata object
  unlink(path, recursive = TRUE)
  snpdata
}

#' Impute missing genotypes for one SNP across all samples
#'
#' @param genotype A vector of integers representing the alleles across all
#'    samples at a given locus.
#'
#' @return The input vector where all missing alleles have been imputed.
#' @keywords internal
#'
impute_one_allele <- function(genotype) {
  checkmate::assert_vector(genotype, min.len = 1L, null.ok = FALSE)
  idx <- which(is.na(genotype))
  for (j in idx) {
    ref         <- length(which(genotype == 0L))
    alt         <- length(which(genotype == 1L))
    maf         <- ifelse(ref < alt, ref / (ref + alt), alt / (ref + alt))
    genotype[j] <- statip::rbern(1L, maf)
  }
  genotype
}
