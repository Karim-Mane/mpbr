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
#'   snpdata <- phase_mixed_genotypes(snpdata)
#'  }
#'
#' @export
phase_mixed_genotypes <- function(snpdata, nsim = 100L) {
  checkmate::assert_class(snpdata, "SNPdata", null.ok = FALSE)
  checkmate::assert_numeric(nsim, lower = 1L, any.missing = FALSE,
                            null.ok = FALSE, len = 1L)
  vcf          <- snpdata[["vcf"]]
  expression   <- "%CHROM\t%POS[\t%AD]\n"
  tmp          <- file.path(dirname(vcf), "tmp")
  system(sprintf("mkdir -p %s", tmp))
  ad           <- file.path(tmp, "AllelicDepth.txt")
  system(sprintf("bcftools query -f'%s' %s > %s", expression, vcf, ad))
  depth        <- data.table::fread(ad, nThread = 4L)
  depth        <- as.matrix(subset(depth, select = -(1L:2L)))
  path         <- file.path(dirname(vcf), "phasing")
  system(sprintf("mkdir -p %s", path))
  correlations <- numeric(length = nsim)
  pb           <- utils::txtProgressBar(min = 0L, max = nsim, initial = 0L,
                                        style = 3L, char = "*")
  for (i in 1L:nsim) {
    tmp_snpdata <- snpdata
    mat         <- apply(tmp_snpdata[["GT"]], 1L, phase_data, depth = depth)
    tmp_snpdata[["Phased"]] <- t(mat)
    saveRDS(t(mat), file.path(path, paste0("sim", i, ".RDS")))
    res_snpdata <- compute_maf(tmp_snpdata,
                               include_het = FALSE,
                               mat_name = "Phased")
    correlations[i] <- stats::cor(res_snpdata[["details"]][["MAF_Phased"]],
                                  res_snpdata[["details"]][["MAF"]])
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
  idx <- which(correlations == max(correlations, na.rm = TRUE))
  snpdata[["Phased"]] <- readRDS(file.path(path,
                                           paste0("sim", idx[[1L]], ".RDS")))
  system(sprintf("rm -rf %s", path))
  snpdata
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

#' Phase the mixed genotypes
#'
#' @param genotype a vector of genotype data
#' @param depth a vector of allelic depth
#'
#' @return a vector of phased genotypes
#' @keywords internal
#'
phase_data <- function(genotype, depth) {
  checkmate::assert_vector(genotype, min.len = 1L, null.ok = FALSE)
  checkmate::assert_vector(depth, min.len = 1L, null.ok = FALSE)
  idx <- as.numeric(which(genotype == 2L))

  for (j in idx) {
    ref <- as.numeric(unlist(strsplit(depth[j], ",", fixed = TRUE))[[1L]])
    alt <- as.numeric(unlist(strsplit(depth[j], ",", fixed = TRUE))[[2L]])
    if (ref == 0L && alt == 0L) {
      ref_count <- sum(genotype == 0L, na.rm = TRUE)
      alt_count <- sum(genotype == 1L, na.rm = TRUE)
      genotype[j] <- phase_ref_alt_0(ref_count, alt_count)
    } else if (ref != 0L && alt != 0L) {
      genotype[j] <- phase_ref_alt_non_0(genotype, ref, alt)
    } else {
      genotype[j] <- phase_ref_or_alt_0(genotype, ref, alt)
    }
  }
  genotype
}
