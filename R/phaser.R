#' Perform phasing of a given mixed genotype
#'
#' @param mad A `data.frame` with three columns that represents the indices of
#'    the mixed alleles, the indices of the SNPs they are found on, and their
#'    allelic depths supporting both reads.
#' @param details A data frame with the SNPs genomic coordinates and their
#'    annotations.
#' @inheritParams phase
#'
#' @returns The input data frame with an additional column named as
#'    `phased_genotype`.
#' @keywords internal
phaser <- function(mad, method, ncores = 4, details = NULL) {
  checkmate::assert_data_frame(mad, min.rows = 1L, ncols = 3L,
                               null.ok = FALSE)
  checkmate::assert_data_frame(details, min.rows = 1L, min.cols = 9L,
                               null.ok = FALSE)
  # ad represents the allelic depth
  mad[["phased_genotype"]] <- switch(
    method,
    major_call = major_call_phaser(mad),
    bi_allelic = bi_allelic_phaser(mad, ncores = ncores),
    least_frequent = least_frequent_phaser(
      mad = mad,
      ncores = ncores,
      details = details
    )
  )
  
  return(mad)
}

#' Phase mixed alleles using the 'major_call' method
#'
#' @inheritParams phaser
#'
#' @inheritSection phaser returns
#' @keywords internal
#'
major_call_phaser <- function(mad, ncores = 4) {
  phase_one_allele <- function(ad) {
    # return NA when the AD is missing
    if (is.na(ad)) {
      return(NA)
    }
    ad <- as.numeric(unlist(strsplit(ad, ",", fixed = TRUE)))
    return(which.max(ad) - 1)
  }
  res <- as.numeric(
    parallel::mclapply(mad[["ad"]], phase_one_allele, mc.cores = ncores)
  )
  return(res)
}

#' Phase mixed alleles using the 'bi_allelic' method
#'
#' @inheritParams phaser
#'
#' @inheritSection phaser returns
#' @keywords internal
#'
bi_allelic_phaser <- function(mad, ncores) {
  # QUESTION: this approach transforms the mixed alleles into either an
  # alternative allele or heterozygous allele regardless of whether the
  # reference allele has more supporting reads.
  # what is the rational behind this?
  phase_one_allele <- function(ad) {
    # return NA when the AD is missing
    if (is.na(ad)) {
      return(NA)
    }
    ad <- as.numeric(unlist(strsplit(ad, ",", fixed = TRUE)))
    read_ratio <- min(ad) / max(ad)
    is_heterozygous <- min(ad) >= 2 & read_ratio > 0.1
    res <- ifelse(is_heterozygous, 2, ifelse(which.max(ad) == 1, 0, 1))
    return(res)
  }
  res <- as.numeric(
    parallel::mclapply(mad[["ad"]], phase_one_allele, mc.cores = ncores)
  )
  return(res)
}

#' Phase mixed alleles using the 'least_frequent' method
#'
#' @inheritParams phaser
#'
#' @inheritSection phaser returns
#' @keywords internal
#'
least_frequent_phaser <- function(mad, ncores, details) {
  # row_index is the row of the corresponding variant. It is derived from:
  # idx_mixed_alleles
  # minor_alle is the code of the minor allele at that variant position. It is
  # derived from the details table obtained after MAF calculation.
  
  phase_one_allele <- function(col, details) {
    # COMBAK: how to deal with cases where the reference and alternate alleles
    # have the same frequencies (MAF_allele = "0=1") ?
    return(details[["MAF_allele"]][col])
  }
  res <- as.numeric(
    parallel::mclapply(
      mad[["col"]], phase_one_allele, details, mc.cores = ncores
    )
  )
  return(res)
}

