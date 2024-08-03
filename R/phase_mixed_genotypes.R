#' Phase mixed genotypes
#'
#' Mixed genotype phasing is performed based on the number of reads supporting
#' each allele at an heterozygous site. The phasing will be run `nsim` times and
#' phased data with the highest correlation between MAF from raw data and MAF
#' from phased data will be retained.
#'
#' @param snpdata a `SNPdata` object
#' @param nsim an integer that represents the number of simulations to be
#'    performed
#' @param genotype The name of the genotype matrix from which the mixed
#'    genotypes will be phased. This can be either the raw genotype matrix
#'    (`GT`) or the imputed genotype matrix (`Imputed`) or any name given to
#'    the target genotype matrix.
#'
#' @return a `SNPdata` object with an additional table named as **Phased**. This
#'    object will contain the phased genotype data
#'
#' @details When both alleles are not supported by any read or the total number
#'    of reads supporting both alleles at a given site is < 5, the genotype will
#'    be phased based on a Bernoulli distribution using the MAF as a parameter.
#'    Similarly, when the total number of reads is > 5 and the number of reads
#'    supporting one of the allele is not 2 times the number of the other,
#'    the genotype is phased using a Bernoulli distribution.
#'
#' @examples
#' \dontrun{
#'   snpdata <- phase(
#'     snpdata,
#'     genotype = "GT"
#'   )
#'  }
#'
#' @export
phase <- function(snpdata, genotype = "GT", nsim = 100L) {
  checkmate::assert_class(snpdata, "SNPdata", null.ok = FALSE)
  checkmate::assert_numeric(nsim, lower = 1L, any.missing = FALSE,
                            null.ok = FALSE, len = 1L)
  checkmate::assert_character(genotype, any.missing = FALSE, len = 1L,
                              null.ok = FALSE)

  # Do not proceed if the specified genotype matrix does not exist
  if (!(genotype %in% names(snpdata))) {
    current_gt_matrices <- names(snpdata)[!(names(snpdata) %in% # nolint: object_usage_linter
                                              c("meta", "details", "vcf"))]
    cli::cli_abort(
      c("x" = "The specified genotype matrix {.code genotype} does not exist",
        "i" = "Current genotype matrices are: {.code {current_gt_matrices}}")
    )
  }

  # Do not proceed if the chosen genotype matrix does not contain any mixed
  # allele.
  if (sum(snpdata[[genotype]] == 2L, na.rm = TRUE) == 0L) {
    cli::cli_abort(
      c("x" = "No mixed genotypes found in the specified genotype matrix: 
        {genotype}",
        "i" = "")
    )
  }

  # extract the allelic depth per sample per SNP. This will be use to choose the
  # most represented allele on a given position for a given sample.
  depth <- extract_allelic_depth(snpdata[["vcf"]])

  # create a temporary directory to store temporary files
  path  <- file.path(tempdir(), "phasing")
  dir.create(path)

  # The correlations vector will store the correlation coefficient between the
  # initial MAF and the MAF after every simulation.
  # The phased data from the simulation where the MAF is the closest to the
  # initial MAF will be chosen. If several of them have the same MAF, the first
  # will be retained.
  correlations <- numeric(length = nsim)
  for (i in cli::cli_progress_along(seq_len(nsim),
                                    "Phasing mixed genotypes")) {
    tmp_snpdata <- snpdata
    depth <- as.matrix(depth)
    mat <- apply(tmp_snpdata[[genotype]], 1L, phase_one_allele, depth)
    tmp_snpdata[["Phased"]] <- t(mat)
    saveRDS(tmp_snpdata[["Phased"]], file.path(path, paste0("sim", i, ".RDS")))
    res_snpdata <- compute_maf(tmp_snpdata,
                               include_het = FALSE,
                               genotype = "Phased")
    correlations[i] <- stats::cor(res_snpdata[["details"]][["MAF_Phased"]],
                                  res_snpdata[["details"]][["MAF"]])
  }
  cli::cli_alert_success("\nMixed genotypes were successfully phased.")
  idx <- which(correlations == max(correlations, na.rm = TRUE))
  snpdata[["Phased"]] <- readRDS(file.path(path,
                                           paste0("sim", idx[[1L]], ".RDS")))
  unlink(path, recursive = TRUE)
  snpdata
}


#' Phase the mixed genotypes for one SNP across all samples
#'
#' @param genotype a vector of genotype data
#' @param depth a vector of allelic depth
#'
#' @return a vector of phased genotypes
#' @keywords internal
#'
phase_one_allele <- function(genotype, depth) {
  checkmate::assert_vector(genotype, min.len = 1L, null.ok = FALSE)
  checkmate::assert_vector(depth, min.len = 1L, null.ok = FALSE)
  depth <- as.character(depth)
  idx <- which(genotype == 2L)

  for (j in idx) {
    ref <- as.numeric(unlist(strsplit(depth[j], ",", fixed = TRUE))[[1L]])
    alt <- as.numeric(unlist(strsplit(depth[j], ",", fixed = TRUE))[[2L]])
    if (any(is.na(c(ref, alt)))) {
      genotype[j] <- phase_missing_ref_or_alt(genotype)
    } else if (ref == 0L && alt == 0L) {
      ref_count   <- sum(genotype == 0L, na.rm = TRUE)
      alt_count   <- sum(genotype == 1L, na.rm = TRUE)
      genotype[j] <- phase_ref_alt_0(ref_count, alt_count)
    } else if (ref != 0L && alt != 0L) {
      genotype[j] <- phase_ref_alt_non_0(genotype, ref, alt)
    } else {
      genotype[j] <- phase_ref_or_alt_0(genotype, ref, alt)
    }
  }
  genotype
}

#' Phase genotype when reference and alternate allele counts = 0
#'
#' @param ref_count the reference allele count
#' @param alt_count the alternate allele count
#'
#' @keywords internal
#'
phase_ref_alt_0 <- function(ref_count, alt_count) {
  if (ref_count < alt_count) {
    genotype <- 0L
  } else if (ref_count > alt_count) {
    genotype <- 1L
  } else {
    genotype <- statip::rbern(1L, ref_count / (ref_count + alt_count))
  }
  genotype
}

#' Phase genotype when both reference and alternate allele counts are > 0
#'
#' @param genotype a vector of integers
#' @param ref the population reference allele count on a specific locus
#' @param alt the population alternate allele count on a specific locus
#'
#' @keywords internal
#'
phase_ref_alt_non_0 <- function(genotype, ref, alt) {
  if (ref + alt >= 5L &&
      (ref >= (2L * alt) || alt >= (2L * ref))) {
    if (ref < alt) {
      res <- 0L
    } else if (ref > alt) {
      res <- 1L
    } else {
      ref_count <- sum(genotype == 0L, na.rm = TRUE)
      alt_count <- sum(genotype == 1L, na.rm = TRUE)
      res       <- phase_ref_alt_0(ref_count, alt_count)
    }
  } else {
    ref_count <- sum(genotype == 0L, na.rm = TRUE)
    alt_count <- sum(genotype == 1L, na.rm = TRUE)
    res       <- phase_ref_alt_0(ref_count, alt_count)
  }
  res
}

#' Phase genotype when either the reference or alternate allele counts = 0
#'
#' @param genotype a vector of integers
#' @param ref the population reference allele count on a specific locus
#' @param alt the population alternate allele count on a specific locus
#'
#' @keywords internal
#'
phase_ref_or_alt_0 <- function(genotype, ref, alt) {
  ref_count <- sum(genotype == 0L, na.rm = TRUE)
  alt_count <- sum(genotype == 1L, na.rm = TRUE)
  if (ref == 0L && alt >= 5L) {
    res     <- 1L
  } else if (ref == 0L && alt < 5L) {
    res     <- statip::rbern(1L, alt_count / (ref_count + alt_count))
  }
  if (alt == 0L && ref >= 5L) {
    res     <- 0L
  } else if (alt == 0L && ref < 5L) {
    res     <- statip::rbern(1L, ref_count / (ref_count + alt_count))
  }
  res
}

#' Phase genotype when either both reference or alternate allelic depth or one
#' of them is missing.
#'
#' @param genotype a vector of integers
#'
#' @keywords internal
#'
phase_missing_ref_or_alt <- function(genotype) {
  ref_count <- sum(genotype == 0L, na.rm = TRUE)
  alt_count <- sum(genotype == 1L, na.rm = TRUE)
  if (ref_count < alt_count) {
    res <- statip::rbern(1L, ref_count / (ref_count + alt_count))
  } else {
    res <- statip::rbern(1L, alt_count / (ref_count + alt_count))
  }
  res
}