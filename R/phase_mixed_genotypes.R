#' Phase mixed genotypes
#'
#' Mixed genotype phasing is performed based on the number of reads supporting
#' each allele at a heterozygous site (allelic depth). A mixed allele will be
#' recoded into either `0` or `1` depending on which one has a highest allelic
#' depth or into `2` when `heterozygous` is set to `TRUE` i.e. to consider
#' mixed alleles as a third allele.
#'
#' @param snpdata a `SNPdata` object
#' @param genotype The name of the genotype matrix from which the mixed
#'    genotypes will be phased. This can be either the raw genotype matrix
#'    (`GT`) or the imputed genotype matrix (`Imputed`) or any name given to
#'    the target genotype matrix.
#' @param ncores A numeric that represents the number of cores to be used during
#'    the phasing process. Default is 4.
#' @param ration A float that represents one of the criteria described below. 
#'    The ratio between the number of reads supporting the two alleles must be
#'    greater than or equal to this cut-off. Default is `0.1`.
#' @param min_ad An integer that represents one of the criteria described below.
#'    The number of reads supporting the reference or alternative alleles must
#'    be greater than or equal to this cut-off. Default is `2`.
#' @param heterozygous A logical that is used to determine whether keep `true`
#'    heterozygous sites (recode them as `2`) or not. Default is `FALSE`. A site
#'    is considered as a true heterozygous when the following two conditions are
#'    met:
#' \enumerate{
#'   \item condition 1: the `ratio` in allelic depth > the specified ratio
#'      cut-off.
#'   \item condition2: the `min_ad` >= the specified minimum allelic depth i.e.
#'      (ratio > 0.1 & min_ad >= 2 by default).
#' }
#' @param phase_equal_ad A character with the name of the phasing approach
#'    used when the allelic depths of the two alleles are equal. The possible
#'    values are:
#' \enumerate{
#'   \item bernoulli: the mixed allele is recoded into `0` or `1` based on a
#'      Bernoulli simulation with a parameter equal to the minor allele
#'      frequency of the SNPs on which the site is found.
#'   \item minor_allele: the mixed allele is transformed into the least
#'      represented allele (i.e. the allele with the minor frequency in the
#'      population).
#' }
#'
#' @return a `SNPdata` object with an additional table named as **Phased** that
#'    contains the phased genotype data.
#'
#' @export
#' 
#' @examples
#' # get SNPdata object
#' snpdata <- get_snpdata(
#'   vcf_file = system.file("extdata", "Input_Data.vcf.gz", package = "mpbr"),
#'   meta_file = system.file("extdata", "SampleMetadata.RDS", package = "mpbr"),
#'   output_dir = tempdir()
#' )
#'
#' # perform mixed genotypes phasing
#' snpdata <- phase(
#'   snpdata = snpdata,
#'   genotype = "GT",
#'   heterozygous = FALSE, # recode all mixed alleles as 0 or 1
#'   ncores = 4,
#'   phase_equal_ad = "bernoulli"
#' )
#'  
phase <- function(snpdata,
                  genotype = "GT",
                  heterozygous = FALSE,
                  min_ad = 2,
                  ratio = 0.1,
                  ncores = 4,
                  phase_equal_ad = "bernoulli") {
  checkmate::assert_class(snpdata, "SNPdata", null.ok = FALSE)
  checkmate::assert_character(genotype, any.missing = FALSE, len = 1L,
                              null.ok = FALSE)
  checkmate::assert_character(phase_equal_ad, any.missing = FALSE, len = 1L,
                              null.ok = FALSE)
  checkmate::assert_logical(
    heterozygous, null.ok = FALSE, any.missing = FALSE, len = 1L
  )
  checkmate::assert_character(phase_equal_ad, any.missing = FALSE, len = 1L,
                              null.ok = FALSE)
  checkmate::assert_choice(
    phase_equal_ad,
    choices = c("bernoulli", "minor_allele"),
    null.ok = FALSE
  )
  checkmate::assert_number(min_ad, na.ok = FALSE, lower = 2, finite = TRUE,
                           null.ok = FALSE)
  checkmate::assert_number(ratio, na.ok = FALSE, finite = TRUE, null.ok = FALSE)
  checkmate::assert_number(ncores, na.ok = FALSE, finite = TRUE,
                           null.ok = FALSE)

  # Do not proceed if the specified genotype matrix does not exist
  check_genotype_matrix(snpdata, genotype)

  # Do not proceed if the chosen genotype matrix does not contain any mixed
  # allele.
  if (sum(snpdata[[genotype]] == 2L, na.rm = TRUE) == 0L) {
    cli::cli_abort(c(
      "x" = "No mixed genotypes found in the specified genotype matrix: \\\
      {genotype}",
      "i" = "The genotype matrix must contain mixed alleles i.e. {.strong 2}.")
    )
  }

  # extract the allelic depth per sample per SNP. This will be use to choose the
  # most represented allele on a given position for a given sample.
  depth <- extract_genotype(snpdata[["vcf"]], which = "AD")
  
  # get the sample ids. This is repeated here to anticipate on the effect of the
  # filtration process. If the phasing is performed after filtration, the
  # remaining samples might not be the same as from the vcf file where the
  # allelic depth is extracted.
  sample_ids <- get_sample_ids(vcf_file = snpdata[["vcf"]])
  
  # The genotype matrix might have undergone some filtration (removal of poor
  # quality samples and SNPs). We filter the allelic depth matrix to only keep
  # the samples and SNPs that are currently in the genotype matrix.
  names(depth) <- c("Chrom", "Pos", sample_ids)
  ad_chrom_pos <- paste(depth[["Chrom"]], depth[["Pos"]], sep = "_")
  gt_chrom_pos <- paste(
    snpdata[["details"]][["Chrom"]],
    snpdata[["details"]][["Pos"]],
    sep = "_"
  )
  row_matches <- ad_chrom_pos %in% gt_chrom_pos
  depth <- depth[row_matches, ]
  col_matches <- which(sample_ids %in% colnames(snpdata[[genotype]]))
  depth <- subset(depth, select = col_matches)
  if (any(dim(depth) != dim(snpdata[[genotype]]))) {
    cli::cli_abort(c(
      x = "The genotype and allelic depth matrices have different dimensions."
    ))
  }
  
  # get indices of the mixed genotypes from allelic depth matrix
  snpdata[["PH"]] <- t(snpdata[[genotype]])
  x <- seq(
    from = 1,
    to = length(snpdata[["PH"]]),
    by = dim(snpdata[["PH"]])[[1L]]
  )
  idx_mixed_alleles <- which(snpdata[["PH"]] == 2L)
  
  # The data has been transposed. Now we have the SNPs in column. We now know
  # that snps on the first column ([line 1 to the end, 1]) belong to first locus
  # Indices of mixed allele which fall within this interval ([1-num_rows]) will
  # be associated to 1 from the `findInterval()` function. That way, we easily
  # pick the MAF of the first SNP and its minor allele to perform the phasing of
  # all mixed alleles on that locus. The same procedure is applied across the
  # rest of the columns.
  transposed_depth <- t(depth[, -c(1:2)])
  intervals <- findInterval(idx_mixed_alleles, x)
  mixed_allele_data <- data.frame(
    idx = idx_mixed_alleles,
    col = intervals,
    ad = transposed_depth[idx_mixed_alleles],
    new_allele = NA
  )
  rm(x)
  
  # identify mixed alleles with no AD (AD=NA) and set them their corresponding
  # genotypes to NA as phasing cannot be performed on them
  if (anyNA(mixed_allele_data[["ad"]])) {
    idx <- mixed_allele_data[["idx"]][is.na(mixed_allele_data[["ad"]])]
    cli::cli_inform(c(
      "!" = "Found {.strong {length(idx)}} mixed alleles with missing allelic depth.",
      "i" = "They will be set to {.strong NA} in the phased genotype matrix."
    ))
    snpdata[["PH"]][idx] <- NA
    mixed_allele_data <- mixed_allele_data[!is.na(mixed_allele_data[["ad"]]), ]
  }
  
  # extract the reference allele count from the allelic depth
  mixed_allele_data[["ref_count"]] <- as.numeric(
    unlist(
      parallel::mclapply(
        mixed_allele_data[["ad"]],
        split_and_extract,
        sep = ",",
        part = 1L,
        mc.cores = ncores
      )
    )
  )
  
  # extract the alternate allele count from the allelic depth
  mixed_allele_data[["alt_count"]] <- as.numeric(
    unlist(
      parallel::mclapply(
        mixed_allele_data[["ad"]],
        split_and_extract,
        sep = ",",
        part = 2L,
        mc.cores = ncores
      )
    )
  )
  
  # get the sum of allelic depths
  mixed_allele_data[["sum"]] <-
    mixed_allele_data[["alt_count"]] + mixed_allele_data[["ref_count"]]
  
  # Determine whether the allelic depths of the reference and alternate alleles
  # are the same. 
  mixed_allele_data[["is_equal"]] <- mixed_allele_data[["ref_count"]] ==
    mixed_allele_data[["alt_count"]]
  
  # get the minimum and maximum AD and 
  mixed_allele_data[["min_ad"]] <- pmin(
    mixed_allele_data[["ref_count"]], mixed_allele_data[["alt_count"]]
  )
  mixed_allele_data[["max_ad"]] <- pmax(
    mixed_allele_data[["ref_count"]], mixed_allele_data[["alt_count"]]
  )
  
  # identify 'true' mixed allele and recode them '2' when 'heterozygous = TRUE'
  if (heterozygous) {
    eff_ratio <- mixed_allele_data[["min_ad"]] / mixed_allele_data[["max_ad"]]
    are_true_mixed <- eff_ratio > ratio &
      mixed_allele_data[["min_ad"]] >= min_ad
    mixed_allele_data[["new_allele"]][are_true_mixed] <- 2
  }
  
  # phase the remaining mixed allele
  if (anyNA(mixed_allele_data[["new_allele"]])) {
    not_same_ad <- is.na(mixed_allele_data[["new_allele"]]) &
      !mixed_allele_data[["is_equal"]]
    if (sum(not_same_ad) > 0) {
      cols <- c("ref_count", "alt_count")
      m <- as.matrix(mixed_allele_data[which(not_same_ad), cols])
      mixed_allele_data[["new_allele"]][which(not_same_ad)] <- max.col(
        replace(m, is.na(m), -Inf), ties.method = "first") - 1
    }
    
    # Now deal with mixed alleles with equal AD
    have_equal_ad <- is.na(mixed_allele_data[["new_allele"]]) &
      mixed_allele_data[["is_equal"]]
    if (sum(have_equal_ad) > 0) {
      idx <- which(have_equal_ad)
      if (phase_equal_ad == "bernoulli") {
        probs <- mixed_allele_data[["min_ad"]][idx] /
          mixed_allele_data[["sum"]][idx]
        mixed_allele_data[["new_allele"]][idx] <- rbinom(
          n = sum(have_equal_ad),
          size = 1,
          prob =  probs
        )
      } else {
        # stop if the minor allele frequencies have not been computed
        if (!("MAF" %in% names(snpdata[["details"]]))) {
          cli::cli_abort(c(
            x = "{.emph minor_allele} method requires minor allele frequency \\\
            to be computed prior to phasing.",
            "!" = "No column named as {.field MAF} found in the \\\
            {.strong details} table of the input {.cls SNPdata} object.",
            i = "Use the {.fn compute_maf} function to calculate the minor \\\
            allele frequency before phasing the data."
          ))
        }
        
        # proceed with phasing if the MAF is already computed
        mixed_allele_data[["new_allele"]][idx] <- as.numeric(
          snpdata[["details"]][["MAF_allele"]][mixed_allele_data[["col"]][idx]]
        )
      }
      
    }
  }
  
  snpdata[["PH"]] <- t(snpdata[[genotype]])
  snpdata[["PH"]][mixed_allele_data[["idx"]]] <- mixed_allele_data[["new_allele"]]
  snpdata[["PH"]] <- t(snpdata[["PH"]])
  
  # the ratio between the AD of the less represented
  # allele and the sum of allelic depths
  # mixed_allele_data[["probs"]] <- unique(
  #   mixed_allele_data[["min_ad"]] / mixed_allele_data[["sum"]]
  # )
  
    # mixed_allele_data[["new_allele"]][are_equal] = rbinom(
    #   n = sum(are_equal),
    #   size = 1,
    #   prob =  mixed_allele_data[["probs"]][are_equal]
    # )
    
  
  
  # when the method is 'least_frequent', check if the minor allele frequency
  # was already computed.
  # Stop the process and ask the user to run 'compute_maf()' on the SNPdata
  # object.
  # if (alternative_method == "minor_allele" &
  #     !("MAF" %in% names(snpdata[["details"]]))) {
  #   cli::cli_abort(c(
  #     x = "{.emph minor_allele} method requires minor allele frequency to \\\
  #     be computed prior to phasing.",
  #     "!" = "No column named as {.field MAF} found in the {.strong details} \\\
  #     table of the input {.cls SNPdata} object.",
  #     i = "Use the {.fn compute_maf} function to calculate the minor allele \\\
  #     frequency before phasing the data."
  #   ))
  # }
  
  # perform the mixed genotypes phasing
  # mixed_allele_data <- phaser(
  #   mad = mixed_allele_data,
  #   method = method,
  #   ratio = ratio,
  #   min_ad = min_ad,
  #   alternative_method = alternative_method,
  #   ncores = ncores,
  #   details = snpdata[["details"]]
  # )

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


phase_genotypes <- function(mixed_allele_data, nsim) {
  checkmate::assert_data_frame(mixed_allele_data, min.rows = 1L, ncols = 3L,
                               null.ok = FALSE)
  # function to split allelic depth on ","
  split_and_extract <- function(x, sep, part) {
    if (is.na(x)) {
      return(NA)
    }
    unlist(strsplit(x, split = sep, fixed = TRUE))[[part]]
  }
  
  # extract the reference allele count from the allelic depth
  mixed_allele_data[["ref_count"]] <- as.numeric(
    unlist(
      parallel::mclapply(
        mixed_allele_data[["ad"]],
        split_and_extract,
        sep = ",",
        part = 1L,
        mc.cores = 4L
      )
    )
  )
  
  # extract the alternate allele count from the allelic depth
  mixed_allele_data[["alt_count"]] <- as.numeric(
    unlist(
      parallel::mclapply(
        mixed_allele_data[["ad"]],
        split_and_extract,
        sep = ",",
        part = 2L,
        mc.cores = 4L
      )
    )
  )
  
  # calculate the difference in allele count between reference and alternate
  # alleles
  mixed_allele_data[["diff"]] <- mixed_allele_data[["ref_count"]] -
    mixed_allele_data[["alt_count"]]
  
  # position where the difference > 1 standard deviation of the difference will
  # be considered as reference. The ones < -1 standard deviation of the
  # difference will be considered as alternate.
  idx <- mixed_allele_data[["diff"]] > sd(mixed_allele_data[["diff"]]) &
    mixed_allele_data[["diff"]]
  
  # send message around the allelic depth in the dataset
  ref_is_major <- round(
    (sum(mixed_allele_data[["diff"]] > 0, na.rm = TRUE) /
       nrow(mixed_allele_data)) * 100,
    digits = 2
  )
  alt_is_major <- round(
    (sum(mixed_allele_data[["diff"]] < 0, na.rm = TRUE) /
       nrow(mixed_allele_data)) * 100,
    digits = 2
  )
  alt_same_as_ref <- round(
    (sum(mixed_allele_data[["diff"]] == 0, na.rm = TRUE) /
       nrow(mixed_allele_data)) * 100,
    digits = 2
  )
  message("Reference allele count is higher in ", ref_is_major,
  "% of the mixed alleles.")
  message("Alternate allele count is higher in ", alt_is_major,
  "% of the mixed alleles.")
  message("Alternate allele count is same as reference count in ",
          alt_same_as_ref,
          "% of the mixed alleles.")
  
  # detect outliers
  y <- mixed_allele_data[["diff"]]
  test_out <- abs(y - median(y, na.rm = TRUE)) / mad(y, na.rm = TRUE)
  # do not plot the below within the function. Return y and test_out and
  # plot at the end of the phasing process
  #plot(y, test_out, ylim=c(0,30), xlim=c(-600,300))
  
  # create a temporary directory to store temporary files
  path  <- file.path(tempdir(), "phasing")
  dir.create(path)
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


