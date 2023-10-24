#' Calculate Weir & Cockerham's Fst
#'
#' @param snpdata `SNPdata` object
#' @param groups a vector of population names. Every sample in the metadata file
#'    is associated to its population of origin. The differentiation index will
#'    be estimated between these groups.
#' @param from the metadata column that contains the sample's population
#'    information.
#'
#' @return a `SNPdata` object with an extra field named as **Fst**. This is a
#'    `list` of data frames with the Fst values for each pair of comparison.
#'    Every data frame contains N rows (where N is the number of loci) and the
#'    following 7 columns:
#'    \enumerate{
#'      \item the chromosome ID
#'      \item the SNPs positions
#'      \item the allele frequency in the first group
#'      \item the allele frequency in the second group
#'      \item the resulting Fst values
#'      \item the p-values associated with the Fst results
#'      \item the p-values corrected for multiple testing using the
#'            Benjamini-Hochberg method
#'    }
#'
#' @examples
#' \dontrun{
#'   snpdata <- calculate_wcFst(
#'     snpdata,
#'     groups = c("Senegal", "Gambia"),
#'     from   = "Country"
#'   )
#' }
#' @export
#'
calculate_wcFst <- function(snpdata, from, groups) { # nolint: object_name_linter
  checkmate::assert_vector(groups, any.missing = FALSE, min.len = 0L,
                           null.ok = TRUE)
  checkmate::assert_character(from, len = 1L, any.missing = FALSE,
                              null.ok = FALSE)
  checkmate::assert_class(snpdata, "SNPdata", null.ok = FALSE)

  if (is.null(groups) && !is.null(from)) {
    groups <- as.character(unique(snpdata[["meta"]][[from]]))
  } else if (all(!is.null(groups) &&
                   !is.null(from) &&
                   (!all(groups %in% unique(snpdata[["meta"]][[from]]))))) {
    stop(sprintf("Not all specified groups belong to the %s column of
                 the metadata table", from))
  }
  system(sprintf("tabix %s", snpdata[["vcf"]]))
  fst <- list()
  for (i in ((seq_len(groups) - 1L))) {
    idx1 <- which(snpdata[["meta"]][[from]] == groups[i])
    idx1 <- paste(idx1, collapse = ",")
    for (j in (i + 1L):length(groups)) {
      idx2 <- which(snpdata[["meta"]][[from]] == groups[j])
      idx2 <- paste(idx2, collapse = ",")
      out  <- file.path(dirname(snpdata[["vcf"]]), "out.wc.fst")
      pout <- file.path(dirname(snpdata[["vcf"]]), "out.wc.fst.pvalues")
      system(sprintf("wcFst --target %s --background %s --file %s \
                     --deltaaf 0 --type GT > %s",
                     idx1, idx2, snpdata[["vcf"]], out))
      system(sprintf("pFst --target %s --background %s --file %s \
                     --deltaaf 0 --type GT > %s",
                     idx1, idx2, snpdata[["vcf"]], pout))
      out        <- data.table::fread(out)
      pout       <- data.table::fread(pout)
      data.table::setkeyv(out,  c("V1", "V2"))
      data.table::setkeyv(pout, c("V1", "V2"))
      tmp        <- out[pout, nomatch = 0L]
      names(tmp) <- c("Chrom", "Pos", "AF_in_target", "AF_in_background",
                      "wcFst", "wcFst_pvalue")
      idx        <- which(tmp[["wcFst"]] < 0L)
      if (length(idx) > 0L) {
        tmp[["wcFst"]][idx] <- 0L
      }
      tmp[["wcFst_Adj_pvalue_BH"]] <- p.adjust(tmp[["wcFst_pvalue"]],
                                               method = "BH")
      fst[[paste0(groups[i], "_vs_", groups[j])]] <- tmp
    }
  }
  snpdata[["Fst"]] <- fst
  snpdata
}

#' Calculate LD R^2 between pairs of loci
#'
#' @param snpdata a `SNPdata` object
#' @param min_r2 the minimum r2 value below which the LD value is not reported
#' @param inter_chrom a `logical` that determines whether to calculate the
#'    inter-chromosomal LD or not. Default is `FALSE`.
#' @param chroms a character vector of chromosome names. If provided, LD will be
#'    calculated only across these chromosomes.
#'
#' @return a `SNPdata` object with an extra field named as **LD**
#' @examples
#' \dontrun{
#'   snpdata <- calculate_ld(snpdata, min_r2=0.2, inter_chrom=FALSE,
#'     chroms=c("Pf3D7_04_v3","Pf3D7_05_v3"))
#'  }
#' @details The output file from LD calculation could be large. In order to
#'     reduce the size of that file, we recommend to specify the list of
#'     chromosomes for which LD should be calculated using the `chroms` option.
#' @export
#'
calculate_ld <- function(snpdata,
                         min_r2      = 0.2,
                         inter_chrom = FALSE,
                         chroms      = NULL) {
  checkmate::assert_class(snpdata, "SNPdata", null.ok = FALSE)
  checkmate::assert_numeric(min_r2, lower = 0L, upper = 1L, finite = TRUE,
                            any.missing = FALSE, null.ok = FALSE)
  checkmate::assert_logical(inter_chrom, any.missing = FALSE, len = 1L,
                            null.ok = FALSE)
  out <- file.path(dirname(snpdata[["vcf"]]), "tmp_ld")
  if (inter_chrom) {
    cat("inter-chromosomal LD will be calculated between sites on",
        paste(chroms, collapse = ","))
    system(sprintf("vcftools --gzvcf %s --out %s --min-r2 %s --interchrom-geno-r2", # nolint: line_length_linter
                   snpdata[["vcf"]], out, min_r2))
  } else {
    system(sprintf("vcftools --gzvcf %s --out %s --min-r2 %s --geno-r2",
                   snpdata[["vcf"]], out, min_r2))
  }
  out <- paste0(out, ".geno.ld")
  system(sprintf("bgzip %s", out))
  out <- paste0(out, ".gz")
  ld  <- data.table::fread(out, nThread = 4L)
  if (!is.null(chroms)) {
    idx1 <- which(ld[["CHR1"]] %in% chroms)
    idx2 <- which(ld[["CHR2"]] %in% chroms)
    idx  <- unique(c(idx1, idx2))
    ld   <- ld[idx, ]
  }

  snpdata[["LD"]] <- ld
  snpdata
}

#' Generate dissimilarity matrix (1-IBS) between all pairs of isolates
#'
#' @param snpdata SNPdata object
#' @param mat_name the name of the genotype table to be used. default="GT"
#'
#' @return a `SNPdata` object with an extra field: `IBS`
#' @examples
#' \dontrun{
#'   snpdata <- calculate_IBS(snpdata, mat_name = "GT")
#'  }
#'
#' @export
#'
calculate_IBS <- function(snpdata, mat_name = "GT") { # nolint: object_name_linter
  if (!(mat_name %in% names(snpdata))) {
    stop("specified genotype matrix does not exist.")
  }
  X  <- t(snpdata[[mat_name]]) # nolint: object_name_linter
  y  <- matrix(NA, nrow = nrow(X), ncol = nrow(X))
  pb <- txtProgressBar(min = 0L, max = nrow(X),
                       initial = 0L, style = 3L, char = "*")
  for (i in seq_along(X)) {
    for (j in seq_along(X)) {
      m       <- X[i, ] - X[j, ]
      y[i, j] <- 1.0 - (length(which(m == 0L)) / ncol(X))
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  y[upper.tri(y)]  <- NA
  colnames(y)      <- rownames(X)
  rownames(y)      <- rownames(X)
  snpdata[["IBS"]] <- y
  snpdata
}


#' Calculate iR index to detect loci with excess of IBD sharing
#'
#' @param snpdata SNPdata object
#' @param mat_name the name of the genotype table to be used. default = "Phased"
#' @param family the name of the column, in the metadata table, to be used to
#'    represent the sample's population
#' @param number_cores the number of cores to be used. default = 4
#'
#' @return a `SNPdata` object with an extra field: `iR`
#' @examples
#' \dontrun{
#'   snpdata <- calculate_iR(
#'     snpdata,
#'     mat_name     = "Phased",
#'     family       = "Location",
#'     number_cores = 4
#'  )
#'  }
#' @export
calculate_iR <- function(snpdata, # nolint: object_name_linter
                         mat_name     = "Phased",
                         family       = "Location",
                         number_cores = 4L) {
  if (!(family %in% names(snpdata[["meta"]]))) {
    stop("No column name ", family, " in the metadata table")
  }
  if (mat_name == "GT" || mat_name == "Phased") {
    cat("Phasing the mixed genotypes\n")
    snpdata <- phase_mixed_genotypes(snpdata, nsim = 10L)
  }
  ped          <- make_ped(snpdata[[mat_name]], snpdata[["meta"]], family)
  ped[["sex"]] <- as.numeric(ped[["sex"]])
  map          <- make_map(snpdata[["details"]])
  ped_map      <- list(ped, map)
  my_geno      <- get_genotypes(ped_map                = ped_map,
                                reference_ped_map      = NULL,
                                maf                    = 0.01,
                                isolate_max_missing    = 0.2,
                                snp_max_missing        = 0.2,
                                chromosomes            = NULL,
                                input_map_distance     = "cM",
                                reference_map_distance = "cM")
  my_param <- isoRelate::getIBDparameters(ped.genotypes = my_geno,
                                          number.cores  = number_cores)
  my_ibd   <- isoRelate::getIBDsegments(ped.genotypes     = my_geno,
                                        parameters        = my_param,
                                        number.cores      = number_cores,
                                        minimum.snps      = 20L,
                                        minimum.length.bp = 50000L,
                                        error             = 0.001)
  my_matrix <- isoRelate::getIBDmatrix(ped.genotypes = my_geno,
                                       ibd.segments  = my_ibd)
  my_iR     <- isoRelate::getIBDiR(ped.genotypes = my_geno, # nolint: object_name_linter
                                   ibd.matrix    = my_matrix,
                                   groups        = NULL)
  pvalues   <- as.numeric(as.character(lapply(my_iR[["log10_pvalue"]],
                                              get_pvalue)))
  my_iR[["log10_pvalue"]]  <- pvalues # nolint: object_name_linter
  names(my_iR)[[8L]]       <- "pvalues" # nolint: object_name_linter
  my_iR[["adj_pvalue_BH"]] <- p.adjust(my_iR[["pvalues"]], method = "BH") # nolint: object_name_linter
  if (!("iR" %in% names(snpdata))) {
    snpdata[["iR"]]        <- list()
  }
  groups    <- unique(snpdata[["meta"]][[family]])
  snpdata[["iR"]][[paste0(groups[[1L]], "_vs_", groups[[2L]])]] <- my_iR
  snpdata
}

get_pvalue <- function(x) {
  10^-x # nolint: implicit_integer_linter
}

make_ped <- function(mat, metadata, family) {
  mat[mat == 1L]  <- 2L
  mat[mat == 0L]  <- 1L
  mat[is.na(mat)] <- 0L
  mat             <- t(mat)

  new_genotype    <- matrix(-9L, nrow = dim(mat)[[1L]],
                            ncol = ((dim(mat)[[2L]]) * 2L))
  j <- 1L
  k <- 1L
  while (j <= ncol(mat)) {
    new_genotype[, k]        <- mat[, j]
    new_genotype[, (k + 1L)] <- mat[, j]
    j                        <- j + 1L
    k                        <- k + 2L
  }

  ped6 <- as.data.frame(cbind(metadata[[family]],
                              metadata[["sample"]],
                              fatherID  = 0L,
                              motherID  = 0L,
                              sex       = metadata[["COI"]],
                              phenotype = -9L),
                        stringsAsFactors = FALSE)
  names(ped6)[1L:2L] <- c("pop", "sample")
  ped                <- as.data.frame(cbind(ped6, data.frame(new_genotype)))
  ped
}

make_map <- function(details) {
  c1      <- details[["Chrom"]]
  c4      <- details[["Pos"]]
  get_xme <- function(x) {
    as.character(unlist(strsplit(x, "_", fixed = TRUE))[[2L]])
  }
  xme     <- as.numeric(as.character(lapply(c1, get_xme)))
  c2      <- paste0(LETTERS[xme], c4)
  c3      <- c4 / 12000L
  map     <- as.data.frame(chrom = c1, post = c2, gd = c3, pos = c4,
                           stringsAsFactors = FALSE)
  map
}

haplotype_to_genotype <- function(haplotypes, moi) {
  genotypes <- matrix(-9L, dim(haplotypes)[[1L]], dim(haplotypes)[[2L]] / 2L)
  for (i in seq_along(dim(haplotypes)[[1L]])) {
    j <- k <- 1L
    while (j <= dim(haplotypes)[[2L]]) {
      if (haplotypes[i, j] == 1L && haplotypes[i, (j + 1L)] == 1L) {
        genotypes[i, k] <- 0L
      }
      if (haplotypes[i, j] == 2L && haplotypes[i, (j + 1L)] == 2L) {
        genotypes[i, k] <- 2L
      }
      if (haplotypes[i, j] == 0L && haplotypes[i, (j + 1L)] == 0L) {
        genotypes[i, k] <- -1L
      }
      if ((moi[i] == 1L) && ((haplotypes[i, j] == 1L && haplotypes[i, (j + 1L)] == 2L) || (haplotypes[i, j] == 2L && haplotypes[i, (j + 1L)] == 1L))) { # nolint: line_length_linter
        genotypes[i, k] <- -1L
      }
      if ((moi[i] == 2L) && ((haplotypes[i, j] == 1L && haplotypes[i, (j + 1L)] == 2L) || (haplotypes[i, j] == 2L && haplotypes[i, (j + 1L)] == 1L))) { # nolint: line_length_linter
        genotypes[i, k] <- 1L
      }
      if (haplotypes[i, j] != 0L && haplotypes[i,j] != 1L && haplotypes[i, j] != 2L) { # nolint: line_length_linter
        genotypes[i, k] <- -1L
      }
      if (haplotypes[i, (j + 1L)] != 0L && haplotypes[i, (j + 1L)] != 1L && haplotypes[i, (j + 1L)] != 2L) { # nolint: line_length_linter
        genotypes[i, k] <- -1L
      }
      j <- j + 2L
      k <- k + 1L
    }
  }
  t(genotypes)
}

calculate_pop_allele_freq <- function(genotypes, moi) {
  pop_allele_freqs <- vector("numeric", dim(genotypes)[[1L]])
  number_isolates  <- dim(genotypes)[[2L]]
  number_snps      <- dim(genotypes)[[1L]]

  for (t in seq_len(number_snps)) {
    A <- 0L # nolint: object_name_linter
    B <- 0L # nolint: object_name_linter
    for (i in seq_len(number_isolates)) {
      if (genotypes[t, i] == 0L && moi[i] == 1L) A <- A + 1L # nolint: object_name_linter
      if (genotypes[t, i] == 0L && moi[i] == 2L) A <- A + 2L # nolint: object_name_linter
      if (genotypes[t, i] == 1L && moi[i] == 2L) { A <- A + 1L; B <- B + 1L } # nolint: object_name_linter
      if (genotypes[t, i] == 2L && moi[i] == 1L) B <- B + 1L # nolint: object_name_linter
      if (genotypes[t, i] == 2L && moi[i] == 2L) B <- B + 1L # nolint: object_name_linter
    }
    if (A + B == 0L) pop_allele_freqs[t] <- -1L # nolint: object_name_linter
    else pop_allele_freqs[t] <- A / (A + B) # nolint: object_name_linter
  }
  pop_allele_freqs
}

calculate_missingness <- function(genotypes) {
  proportion_missing <- vector("numeric", dim(genotypes)[[2L]])
  number_snps        <- dim(genotypes)[[2L]]
  number_isolates    <- dim(genotypes)[[1L]]
  number_snps_1      <- dim(genotypes)[[1L]]

  for (i in seq_len(number_snps)) {
    number_missing <- 0L
    for (j in seq_len(number_isolates)) {
      if (genotypes[j, i] == -1L) number_missing <- number_missing + 1L
    }
    proportion_missing[i] <- number_missing / number_snps_1
  }
  proportion_missing
}

check_input_ped_map_files <- function(input_map, input_ped,
                                      input_map_distance) {
  checkmate::assert_character(input_map_distance, len = 1L, null.ok = FALSE)
  checkmate::assert_choice(input_map_distance, choices = c("M", "cM"))
  stopifnot("'ped_map' has incorrect format - genetic-map positions in MAP file # nolint: consecutive_assertion_linter
            are non-numeric" = is.numeric(input_map[, "pos_M"]))
  stopifnot("'ped_map' has incorrect format - base-pair positions in MAP file
            are non-numeric" = is.numeric(input_map[, "pos_bp"]))

  if (is.factor(input_map[, "chr"])) {
    input_map[, "chr"] <- as.character(input_map[, "chr"])
  }
  if (is.factor(input_map[, "snp_id"])) {
    input_map[, "snp_id"] <- as.character(input_map[, "snp_id"])
  }
  if (is.factor(input_ped[, 1L])) {
    input_ped[, 1L] <- as.character(input_ped[, 1L])
  }
  if (is.factor(input_ped[, 2L])) {
    input_ped[, 2L] <- as.character(input_ped[, 2L])
  }
  if (input_map_distance == "cM") {
    input_map[, "pos_M"] <- input_map[, "pos_M"] / 100L
  }
  list(
    input_map = input_map,
    input_ped = input_ped
  )
}

check_ref_ped_map_files <- function(reference_ped_map, reference_map_distance) {
  checkmate::assert_list(reference_ped_map, len = 2L, null.ok = TRUE)
  checkmate::assert_character(reference_map_distance, len = 1L, null.ok = FALSE)
  checkmate::assert_choice(reference_map_distance, choices = c("M", "cM"))

  reference_ped <- reference_map <- NULL
  if (!is.null(reference_ped_map)) {
    reference_ped <- reference_ped_map[[1L]]
    reference_map <- reference_ped_map[[2L]]
    stopifnot("'reference_ped_map' has incorrect format - PED is not # nolint: consecutive_assertion_linter
              a data.frame" = is.data.frame(reference_ped))
    stopifnot("'reference_ped_map' has incorrect format - MAP is not # nolint: consecutive_assertion_linter
              a data.frame" = is.data.frame(reference_map))
    # check the PED and MAP files have the same number of SNPs
    stopifnot("'reference_ped_map' has incorrect format - PED and MAP must have # nolint: consecutive_assertion_linter
              the same number of SNPs" =
                ncol(reference_ped) == (2L * nrow(reference_map) + 6L))

    # check the MAP file has 4 columns
    stopifnot("'reference_ped_map' has incorrect format - MAP must have
              four columns" = ncol(reference_map) == 4L)
    colnames(reference_map) <- c("chr", "snp_id", "pos_M", "pos_bp")
    stopifnot("'reference_ped_map' has incorrect format - genetic-map positions # nolint: consecutive_assertion_linter
              in reference MAP file are non-numeric" =
                is.numeric(reference_map[, "pos_M"]))
    stopifnot("'reference_ped_map' has incorrect format - base-pair positions
              in reference MAP file are non-numeric" =
                is.numeric(reference_map[, "pos_bp"]))

    # check MAP for factors
    if (is.factor(reference_map[, "chr"])) {
      reference_map[, "chr"] <- as.character(reference_map[, "chr"])
    }
    if (is.factor(reference_map[, "snp_id"])) {
      reference_map[, "snp_id"] <- as.character(reference_map[, "snp_id"])
    }

    # check PED for factors
    if (is.factor(reference_ped[, 1L])) {
      reference_ped[, 1L] <- as.character(reference_ped[, 1L])
    }
    if (is.factor(reference_ped[, 2L])) {
      reference_ped[, 2L] <- as.character(reference_ped[, 2L])
    }
    if (reference_map_distance == "cM") {
      reference_map[, "pos_M"] <- reference_map[, "pos_M"] / 100L
    }
  }

  list(
    reference_ped = reference_ped,
    reference_map = reference_map
  )
}

check_chromosomes <- function(chromosomes, input_map,
                              reference_ped_map, reference_map_distance) {
  checkmate::assert_vector(chromosomes, null.ok = TRUE)
  # check chromosomes
  stopifnot("'chromosomes' has incorrect format - must be a vector" =
              is.vector(chromosomes))
  if (!is.null(chromosomes)) {
    if (!all(chromosomes %in% input_map[, "chr"])) {
      stop(sprintf("chromosome %s not in 'ped_map'\n",
                   paste0(chromosomes[!(chromosomes %in% input_map[, "chr"])])))
    }

    if (!is.null(reference_ped_map)) {
      res           <- check_ref_ped_map_files(reference_ped_map,
                                               reference_map_distance)
      reference_map <- res[["reference_map"]]
      if (!all(chromosomes %in% reference_map[, "chr"])) {
        stop(sprintf("chromosome %s not in 'reference_ped_map'\n",
                     paste0(chromosomes[!(chromosomes %in% reference_map[, "chr"])]))) # nolint: line_length_linter 
      }
    }
  } else {
    chromosomes <- unique(as.character(input_map[, "chr"]))
  }
  chromosomes
}

#' Title
#'
#' @param ped_map 
#' @param reference_ped_map 
#' @param maf 
#' @param isolate_max_missing 
#' @param snp_max_missing 
#' @param chromosomes 
#' @param input_map_distance 
#' @param reference_map_distance 
#'
#' @return
#' @export
#'
#' @examples
get_genotypes <- function(ped_map,
                          reference_ped_map      = NULL,
                          maf                    = 0.01,
                          isolate_max_missing    = 0.1,
                          snp_max_missing        = 0.1,
                          chromosomes            = NULL,
                          input_map_distance     = "cM",
                          reference_map_distance = "cM") {
  checkmate::assert_numeric(maf, lower = 0.0, upper = 1.0, finite = TRUE,
                            null.ok = FALSE)
  checkmate::assert_list(ped_map, len = 2L, null.ok = FALSE)
  checkmate::assert_data_frame(ped_map[[1L]],
                               ncols = (2L * nrow(ped_map[[2L]]) + 6L))
  checkmate::assert_data_frame(ped_map[[2L]], ncols = 4L)
  checkmate::assert_numeric(isolate_max_missing, lower = 0.0, upper = 1.0,
                            finite = TRUE, null.ok = FALSE)
  checkmate::assert_numeric(snp_max_missing, lower = 0.0, upper = 1.0,
                            finite = TRUE, null.ok = FALSE)
  checkmate::assert_character(input_map_distance, len = 1L, null.ok = FALSE)
  checkmate::assert_character(reference_map_distance, len = 1L, null.ok = FALSE)
  checkmate::assert_choice(input_map_distance, choices = c("M", "cM"),
                           null.ok = FALSE)
  checkmate::assert_choice(reference_map_distance, choices = c("M", "cM"),
                           null.ok = FALSE)
  checkmate::assert_list(reference_ped_map, null.ok = TRUE, len = 2L,
                         names = c("PED", "MAP"))

  # check the input ped and map files
  input_ped <- as.data.frame(ped_map[[1L]])
  input_map <- as.data.frame(ped_map[[2L]])
  colnames(input_map) <- c("chr", "snp_id", "pos_M", "pos_bp")
  check_res <- check_input_ped_map_files(input_map, input_ped,
                                         input_map_distance)
  input_ped <- check_res[["input_ped"]]
  input_map <- check_res[["input_map"]]

  # check chromosomes
  chromosomes <- check_chromosomes(chromosomes, input_map,
                                   reference_ped_map, reference_map_distance)

  # create new isolate IDs from PED FIDs and IIDs
  isolate_names <- paste(input_ped[, 1L], input_ped[, 2L], sep = "/")
  if (anyDuplicated(isolate_names) > 0L) {
    stop("\nDuplicate sample IDs found.")
  }

  ## ---
  ## begin data filtering
  ## ---

  # merge input data with reference data
  if (!is.null(reference_ped_map)) {
    # check the reference ped and map files
    check_res <- check_ref_ped_map_files(reference_ped_map,
                                         reference_map_distance)
    reference_ped <- check_res[["reference_ped"]]
    reference_map <- check_res[["reference_map"]]

    input_map_v1      <- cbind(seq_along(input_map), input_map)
    reference_map_v1  <- cbind(seq_along(reference_map), reference_map)
    input_map_v1      <- merge(input_map_v1, reference_map_v1,
                               by.x = "snp_id", by.y = "snp_id")
    if (nrow(input_map_v1) == 0L) {
      stop("\nNo SNPs remaining after merging 'ped_map' and 'reference_ped_map'") # nolint: line_length_linter
    }

    ## ---
    ## need an example to correct the following
    input_map_v1      <- input_map_v1[order(input_map_v1[, "1:nrow(input_map)"]), ] # nolint: line_length_linter
    ## ---

    if (!is.null(chromosomes)) {
      input_map_v1 <- input_map_v1[input_map_v1[, "chr.x"] %in% chromosomes, ]
    }
    stopifnot("No SNPs remaining after merging 'ped_map' and 'reference_ped_map'
              for selected chromosomes" = nrow(input_map_v1) != 0L)

    input_ped_columns <- c(1L:6L, 2L * input_map_v1[, "1:nrow(input_map)"] + 5L,
                           2L * input_map_v1[, "1:nrow(input_map)"] + 6L)
    input_ped_columns <- input_ped_columns[sort(input_ped_columns,
                                                na.last = TRUE)]
    input_ped_v1      <- input_ped[, input_ped_columns]
    reference_ped_columns  <- c(1L:6L, 2L * input_map_v1[, "1:nrow(reference.map)"] + 5L, # nolint: line_length_linter
                                2L * input_map_v1[, "1:nrow(reference.map)"] + 6L) # nolint: line_length_linter
    reference_ped_columns  <- reference_ped_columns[sort(reference_ped_columns,
                                                         na.last = TRUE)]
    reference_ped_v1       <- reference_ped[, reference_ped_columns]
    input_map_v2           <- input_map_v1[, c("chr.x", "snp_id", "pos_M.x",
                                               "pos_bp.x")]
    colnames(input_map_v2) <- c("chr", "snp_id", "pos_M", "pos_bp")
  } else {
    if (!is.null(chromosomes)) {
      input_map_v1 <- cbind(seq_along(input_map), input_map)
      input_map_v1 <- input_map_v1[input_map_v1[, "chr"] %in% chromosomes, ]
      if (nrow(input_map_v1) == 0L) {
        stop("No SNPs remaining after subsetting 'ped_map' by selected
              chromosomes")
      }
      input_ped_columns <- c(1L:6L, 2L * input_map_v1[, "1:nrow(input_map)"] + 5L, # nolint: line_length_linter
                             2L * input_map_v1[, "1:nrow(input_map)"] + 6L)
      input_ped_columns <- input_ped_columns[sort(input_ped_columns,
                                                  na.last = TRUE)]
      input_ped_v1      <- input_ped[, input_ped_columns]
      input_map_v2      <- input_map_v1[, c("chr", "snp_id", "pos_M", "pos_bp")]
    } else {
      input_map_v2 <- input_map
      input_ped_v1 <- input_ped
    }
  }

  # call genotypes
  input_matrix        <- as.matrix(input_ped_v1[, 7L:ncol(input_ped_v1)])
  input_genders       <- input_ped_v1[, 5L]
  input_genotypes_v0  <- cbind(input_map_v2,
                               haplotype_to_genotype(input_matrix,
                                                     input_genders))
  if (!is.null(reference_ped_map)) {
      reference_matrix       <- as.matrix(reference_ped_v1[, 7L:ncol(reference_ped_v1)]) # nolint: line_length_linter
      reference_genders      <- reference_ped_v1[, 5L]
      reference_genotypes_v0 <- cbind(input_map_v2,
                                      haplotype_to_genotype(reference_matrix,
                                                            reference_genders))
  }

  # calculate allele frequencies form reference data
  if (is.null(reference_ped_map)) {
    pop_allele_freq    <- calculate_pop_allele_freq(
      as.matrix(input_genotypes_v0[, 5L:ncol(input_genotypes_v0)]),
      input_ped_v1[, 5L]
    )
    input_genotypes_v1 <- cbind(input_genotypes_v0[, 1L:4L],
                                pop_allele_freq,
                                input_genotypes_v0[, 5L:ncol(input_genotypes_v0)]) # nolint: line_length_linter
  } else {
    pop_allele_freq    <- calculate_pop_allele_freq(
      as.matrix(reference_genotypes_v0[, 5L:ncol(reference_genotypes_v0)]),
      reference_ped_v1[, 5L]
    )
    input_genotypes_v1 <- cbind(input_genotypes_v0[, 1L:4L],
                                pop_allele_freq,
                                input_genotypes_v0[, c(5L:ncol(input_genotypes_v0))]) # nolint: line_length_linter
  }
  colnames(input_genotypes_v1) <- c("chr", "snp_id", "pos_M", "pos_bp", "freq",
                                    isolate_names)
  cat(paste0("\nBegin filtering of", length(isolate_names), "isolates and ",
             nrow(input_genotypes_v1), "SNPs...\n"))

  # remove SNPs with low population MAF
  # removing SNPs with AF>0.99 and AF<0.01
  input_genotypes_v2 <- subset(input_genotypes_v1,
                               pop_allele_freq <= (1.0 - maf) & pop_allele_freq >= maf) # nolint: line_length_linter
  stopifnot("0 SNPs remain after MAF removal" = nrow(input_genotypes_v2) > 0L)
  cat(paste0(nrow(input_genotypes_v2),
             "SNPs remain after MAF removal...\n"))

  # remove snps with high missingness
  snp_missingness    <- calculate_missingness(
    as.matrix(t(input_genotypes_v2[, 6L:ncol(input_genotypes_v2)]))
  )
  input_genotypes_v3 <- input_genotypes_v2[snp_missingness <= snp_max_missing, ]
  stopifnot("0 SNPs remain after missingness removal" =
              nrow(input_genotypes_v3) > 0L)
  cat(paste0(nrow(input_genotypes_v3),
             "SNPs remain after missingness removal...\n"))

  # remove samples with high missingness
  isolate_missingness <- round(
    calculate_missingness(
      as.matrix(input_genotypes_v3[, 6L:ncol(input_genotypes_v3)])
    ), digits = 3L
  )
  if (length(isolate_names[isolate_missingness > isolate_max_missing]) > 0L) {
    my_remove <- isolate_names[isolate_missingness > isolate_max_missing]
    warning("isolates removed due to genotype missingness: ",
            glue::glue_collapse(my_remove, sep = ", "), call. = FALSE)
    sample_keep        <- input_ped_v1[isolate_missingness <= isolate_max_missing, 1L:6L] # nolint: line_length_linter
    input_genotypes_v4 <- input_genotypes_v3[, c(1L:5L, which(isolate_missingness <= isolate_max_missing) + 5L)] # nolint: line_length_linter
    if (nrow(sample_keep) < 1L) {
      stop(
           sprintf("All isolates removed with missingness > %s%s. No isolates remaining.", # nolint: line_length_linter
                   isolate_max_missing * 100.0, "%"))
    }
  } else {
    sample_keep        <- input_ped_v1[, 1L:6L]
    input_genotypes_v4 <- input_genotypes_v3
  }
  colnames(sample_keep) <- c("fid", "iid", "pid", "mid", "moi", "aff")
  if ((ncol(input_genotypes_v4) - 5L) == 0L) {
    stop("0 samples remain after missingness removal")
  }
  cat(paste0(ncol(input_genotypes_v4) - 5L,
             "isolates remain after missingness removal...\n"))

  return_genotypes        <- list(
    pedigree  = sample_keep,
    genotypes = input_genotypes_v4
  )
  return_genotypes
}

#' Estimate the relatedness between pairs of isolates
#' @param snpdata SNPdata object
#' @param mat_name the name of the genotype table to be used. default="GT"
#' @param from the name of the column, in the metadata table, to be used to represent the sample's population
#' @param sweepRegions a data frame with the genomic coordinates of the regions of the genome to be discarded. This should contain the following 3 columns:
#' \enumerate{
#' \item Chrom: the chromosome ID
#' \item Start: the start position of the region on the chromosome
#' \item End: the end position of the region on the chromosome
#' }
#' @param groups a vector of character. If specified, relatedness will be generated between these groups
#' @return SNPdata object with an extra field: relatedness. This will contain the relatedness data frame of 3 columns and its correspondent matrix
#' @examples
#' \dontrun{
#'   calculate_relatedness(snpdata, mat_name="GT", family="Location", 
#'   sweepRegions=NULL, groups=c("Chogen","DongoroBa"))
#'  }
#' @details The relatedness calculation is based on the model developed by Aimee R. Taylor and co-authors. https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009101
#' @export
calculate_relatedness = function(snpdata, mat_name="Imputed", from="Location", sweepRegions=NULL, groups=NULL){
    # sourceCpp("src/hmmloglikelihood.cpp")
    details = snpdata$details
    metadata = snpdata$meta
    if(mat_name=="GT" | mat_name=="Phased"){
        cat("Imputing the missing genotypes\n")
        snpdata = impute_missing_genotypes(snpdata, genotype=mat_name, nsim=10)
    }
    if((mat_name=="Imputed") & (!("Imputed" %in% names(snpdata)))){
        mat_name="GT"
        cat("Imputing the missing genotypes\n")
        snpdata = impute_missing_genotypes(snpdata, genotype=mat_name, nsim=10)
    }

    # mat_name = "Imputed"
    mat = snpdata[["Imputed"]]
    if(!is.null(groups) & all(groups %in% unique(metadata[[from]]))){
        pops = groups
    }else{
        pops = unique(metadata[[from]])
    }
    cat("creating the genotype data\n")
    genotypes = constructGenotypeFile(pops, details, mat, metadata, from, sweepRegions)
    sites = names(genotypes)
    dir = dirname(snpdata$vcf)
    ibd = NULL
    for(ii in 1:length(sites)){
        site1 = sites[ii]
        for(jj in ii:length(sites)){
            site2 = sites[jj]
            cat("calculating the relatedness between",site1,"and",site2,"\n")
            ibd = rbind(ibd, gen_mles(genotypes, site1, site2, f=0.3, dir))
        }
    }
    ibd = as.data.table(ibd)
    names(ibd) = c('iid1','iid2','k','relatedness')
    ibd = subset(ibd, select = -3)
    ibd$relatedness = as.numeric(ibd$relatedness)
    cat("creating the relatedness matrix\n")
    rmatrix = createRelatednessMatrix(ibd, metadata, pops, from)
    snpdata$relatedness = list()
    snpdata$relatedness[["df"]] = ibd
    snpdata$relatedness[["matrix"]] = rmatrix
    system(sprintf("rm -rf %s", paste(dir,"/ibd")))
    snpdata
}

createRelatednessMatrix = function(ibd, metadata, pops, from){
    idx = NULL
    for(pop in pops){
        idx = c(idx, which(metadata[[from]]==pop))
    }
    idx = unique(idx)
    metadata = metadata[idx,]
    modifier = function(x){gsub('-','.',x)}
    metadata$sample = as.character(lapply(metadata$sample,modifier))
    relatednessMatrix = matrix(NA,nrow = length(metadata$sample), ncol = length(metadata$sample))
    rownames(relatednessMatrix) = metadata$sample
    colnames(relatednessMatrix) = metadata$sample
    surLigne = unique(ibd$iid1)
    for(i in 1:length(surLigne)){
        # print(paste0('i=',i))
        ligne = surLigne[i]
        l = match(ligne,rownames(relatednessMatrix))
        target = ibd[which(ibd$iid1==ligne),]
        t = unique(target$iid2)
        for(j in 1:length(t)){
            tt = target[which(target$iid2==t[j]),]
            k = match(t[j],colnames(relatednessMatrix))
            if(nrow(tt) == 0)
                relatednessMatrix[l,k] = 0
            else if(nrow(tt)==1)
                relatednessMatrix[l,k] = round(as.numeric(tt$relatedness), digits = 5)
            else if(nrow(tt)>1 & length(unique(tt$iid1))==1 & length(unique(tt$iid2))==1)
                relatednessMatrix[l,k] = round(as.numeric(tt$relatedness[1]), digits = 5)
            else if(nrow(tt)>1 & length(unique(tt$iid1))>1 | length(unique(tt$iid2))>1)
                relatednessMatrix[l,k] = mean(round(as.numeric(tt$relatedness), digits = 5), na.rm = TRUE)
        }
    }
    relatednessMatrix
}

constructGenotypeFile = function(pops, details, mat, metadata, from, sweepRegions){
    if(!is.null(sweepRegions)){
        selectiveRegions = fread(sweepRegions)
        rtd = NULL
        for(j in 1:nrow(selectiveRegions))
            rtd = c(rtd, which(details$Chrom==selectiveRegions$Chrom[j] & (details$Pos>=selectiveRegions$Start[j] & details$Pos<=selectiveRegions$End[j])))
        rtd = unique(rtd)
        details = details[-rtd,]
        mat = mat[-rtd,]
    }
    res = list()
    # pops = unique(metadata[[from]])
    for(pop in pops){
        idx = which(metadata[[from]]==pop)
        s = metadata$sample[idx]
        m = match(s, colnames(mat))
        X=mat[,m]
        chroms=details$Chrom; pos=details$Pos; samps=metadata$sample[idx]
        L = list(X, chroms, pos, samps)
        names(L)=c('X','chroms','pos','samps')
        res[[pop]] = L
    }
    res
}

gen_mles = function(res, site1, site2, f=0.3, dir){
    nproc=700
    epsilon=0.001
    nboot=100
    Ps=c(0.025,0.975)
    outFilePrefix = paste0(site1,'_',site2)
    countryA = res[[site1]]
    countryB = res[[site2]]
    outputDir = paste0(dir,'/ibd')
    system(sprintf("mkdir -p %s", outputDir))
    if (site1==site2){
        # within country comparison
        L=run_country(countryA,f)
        data_set=L$data_set
        individual_names=L$individual_names
        nindividuals=L$nindividuals
        chrom=L$chrom
        pos=L$pos
        name_combinations = matrix(nrow = nindividuals*(nindividuals-1)/2, ncol = 2)
        Y=matrix(data=NA,nrow =length(name_combinations),ncol = 4 )
        k=0
        for ( i in 1 : (nindividuals-1)){
            j=(1+k):(k+nindividuals-i)
            name_combinations[j,1]=rep(individual_names[i],each=length(j))
            name_combinations[j,2]=individual_names[(i+1):nindividuals]
            k=k+length(j)      # within country comparison
        }
    }else{
        # within country comparison
        LA=run_country(countryA,f)
        individual_names_A=LA$individual_names
        nindividuals_A=LA$nindividuals
        # # within country comparison
        LB=run_country(countryB,f)
        individual_names_B=LB$individual_names
        nindividuals_B=LB$nindividuals
        c_offset = dim(LA$data_set)[2]
        L=run_2country(countryA,countryB,f)
        data_set=L$data_set
        individual_names=L$individual_names
        nindividuals=L$nindividuals
        chrom=L$chrom
        pos=L$pos
        name_combinations <- matrix(nrow = nindividuals_A*nindividuals_B, ncol = 2)
        Y=matrix(data=NA,nrow =length(name_combinations),ncol = 4 )
        name_combinations[,1]=rep(individual_names_A,each=nindividuals_B)
        name_combinations[,2]=rep(individual_names_B,nindividuals_A)
    }
    X=as.matrix(data_set)
    data_set$fs = rowMeans(data_set, na.rm = TRUE) # Calculate frequencies
    data_set$pos =pos
    data_set$chrom=chrom
    data_set$dt <- c(diff(data_set$pos), Inf)
    pos_change_chrom <- 1 + which(diff(data_set$chrom) != 0) # find places where chromosome changes
    data_set$dt[pos_change_chrom-1] <- Inf
    # note NA result is undocumented - could change
    a0=Rfast::rowCountValues(X, rep(0,dim(X)[1]))  #getting the count of 0 on each row (SNPs)
    a1=Rfast::rowCountValues(X, rep(1,dim(X)[1]))  #getting the count of 1 on each row (SNPs)
    a2=Rfast::rowCountValues(X, rep(2,dim(X)[1]))  #getting the count of 2 on each row (SNPs)
    ana=rowSums(is.na(X))
    frequencies=cbind(a0,a1,a2)/(dim(X)[2]-ana)
    if (all.equal(rowSums(frequencies),rep(1,dim(X)[1]))!=TRUE){
        cat(paste0("frequency ERROR"))
        return()
    }
    # if (iloop < nproc){
    #     N=floor(dim(name_combinations)[1]/nproc)
    #     starty=(iloop-1)*N+1
    #     endy=iloop*N
    # }else{
    #     N=dim(name_combinations)[1]-(nproc-1)*floor(dim(name_combinations)[1]/nproc)
    #     starty = 1+(nproc-1)*floor(dim(name_combinations)[1]/nproc)
    #     endy=dim(name_combinations)[1]
    # }
    starty = 1
    endy = dim(name_combinations)[1]
    X=matrix(data=NA,nrow =endy-starty+1,ncol = 4 )
    for (icombination in starty:endy){
        #cat(paste0("icombination=",icombination,"\n"))
        individual1 <- name_combinations[icombination,1]
        individual2 <- name_combinations[icombination,2]
        if (site1==site2){
            # Indices of pair
            i1 = which(individual1 == names(data_set))  #index of ind1 on the data frame
            i2 = which(individual2 == names(data_set)) #index of ind2 on the data frame
        }else{
            i1 = which(individual1 == individual_names_A) #index of ind1 on the data frame
            i2 = which(individual2 == individual_names_B) + c_offset #index of ind2 on the data frame
        }
        # Extract data
        subdata <- cbind(data_set[,c("fs","dt")],data_set[,c(i1,i2)])
        names(subdata) <- c("fs","dt","Yi","Yj") # note fs not used
        krhat_hmm <- compute_rhat_hmm(frequencies, subdata$dt, cbind(subdata$Yi, subdata$Yj), epsilon)
        X[icombination-starty+1,1]=individual1
        X[icombination-starty+1,2]=individual2
        X[icombination-starty+1,3:4]=krhat_hmm

    }
    saveRDS(X, file = paste0(outputDir,'/',outFilePrefix,"_",starty,"_",endy,"_",gsub("\\.","p",sprintf("%.4f",f)),".RDS"))
    X
}

compute_rhat_hmm = function(frequencies, distances, Ys, epsilon){
    ll <- function(k, r) loglikelihood(k, r, Ys, frequencies, distances, epsilon, rho = 7.4 * 10^(-7)) #loglikelihood_cpp(k, r, Ys, frequencies, distances, epsilon, rho = 7.4 * 10^(-7))
    optimization <- optim(par = c(50, 0.5), fn = function(x) - ll(x[1], x[2]))
    rhat <- optimization$par
    return(rhat)
}

## Mechanism to generate Ys given fs, distances, k, r, epsilon
simulate_Ys_hmm <- function(frequencies, distances, k, r, epsilon){
    Ys <- simulate_data(frequencies, distances, k = k, r = r, epsilon, rho = 7.4 * 10^(-7))
    return(Ys)
}


run_country=function(countryA,f){
    matchy=function(x){regmatches(x, regexec('_(.*?)\\_', x))[[1]][2]}
    L=countryA  #readRDS(countryA)
    q=data.frame(L$X)
    q=lapply(q,function(x) as.integer(x))
    Q=matrix(unlist(q), ncol = length(q[[1]]), byrow = TRUE)
    maf=colSums(Q==1,na.rm =TRUE)/colSums(Q<=1,na.rm=TRUE) #colSums(!is.na(Q))
    i=which(colSums(Q==0,na.rm =TRUE)<colSums(Q==1,na.rm =TRUE))
    maf[i]=colSums(Q[,i]==0,na.rm =TRUE)/colSums(Q[,i]<=1,na.rm=TRUE) #colSums(!is.na(Q[,i]))
    j=which(maf<=f)
    data_set = data.frame(q)[-j,]
    # Create indices for pairwise comparisons
    individual_names <- names(q)
    nindividuals <- length(individual_names)
    chrom=as.integer(unlist(lapply(L$chroms[-j],matchy)))
    pos=L$pos[-j]
    return(list(data_set=data_set,j=j,
                individual_names=individual_names,nindividuals=nindividuals,
                chrom=chrom,pos=pos))
}

run_2country=function(countryA,countryB,f){
    matchy=function(x){regmatches(x, regexec('_(.*?)\\_', x))[[1]][2]}
    Z=countryA  #readRDS(countryA)
    M=countryB  #readRDS(countryB)
    X=cbind(Z$X,M$X)
    samps=c(Z$samps,M$samps)
    chroms=c(Z$chroms)
    pos=c(Z$pos)#
    q=data.frame(X)
    q=lapply(q,function(x) as.integer(x))
    Q=matrix(unlist(q), ncol = length(q[[1]]), byrow = TRUE)
    maf=colSums(Q==1,na.rm =TRUE)/colSums(Q<=1,na.rm=TRUE) #colSums(!is.na(Q))
    i=which(colSums(Q==0,na.rm =TRUE)<colSums(Q==1,na.rm =TRUE))
    maf[i]=colSums(Q[,i]==0,na.rm =TRUE)/colSums(Q[,i]<=1,na.rm=TRUE) #colSums(!is.na(Q[,i]))
    j=which(maf<=f)
    data_set = data.frame(q)[-j,]
    # Create indices for pairwise comparisons
    individual_names <- names(q)
    nindividuals <- length(individual_names)
    chrom=as.integer(unlist(lapply(chroms[-j],matchy)))
    pos=pos[-j]
    return(list(data_set=data_set,j=j,
                individual_names=individual_names,nindividuals=nindividuals,
                chrom=chrom,pos=pos))
}

loglikelihood = function(k, r, Ys, f, gendist, epsilon, rho = 7.4 * 10^(-7)){
    loglikelihood_value = 0
    if (r < 0 | r > 1 | k < 0){
        return(log(0))
    }
    current_predictive = numeric(length = 2)
    current_predictive[1] = 1 - r
    current_predictive[2] = r
    current_filter = numeric(length = 2)
    ndata = nrow(Ys)
    maxnstates = ncol(f)
    for (idata in 1:ndata){
        nstates = 1
        # print(paste0("idata=",idata))
        while((nstates <= maxnstates) && (f[idata,nstates] > 1e-20)){
            nstates = nstates+1
        }
        lk0 = 0
        incr = 0
        # print(paste0("nstates=",nstates))
        # print(Ys[idata,])
        if(nstates>maxnstates) nstates=maxnstates
        for (g in 1:nstates){
            for (gprime in 1:nstates){
                if(gprime<=nstates){
                    incr = f[idata, g] * f[idata, gprime]
                    if (Ys[idata,1] == g){
                        incr = incr * (1 - (nstates - 1) * epsilon)
                    } else {
                        incr = incr*epsilon
                    }
                    if (Ys[idata,2] == gprime){
                        incr = incr * (1 - (nstates - 1) * epsilon)
                    } else {
                        incr = incr*epsilon
                    }
                    lk0 = lk0 + incr
                }
            }
        }
        lk1 = 0
        incr = 0

        for (g in 1:nstates){
          incr = f[idata, g]
          if (Ys[idata,1] == g){
            incr = incr * (1 - (nstates - 1) * epsilon)
          } else {
            incr = incr*epsilon
          }
          if (Ys[idata,2] == g){
            incr = incr * (1 - (nstates - 1) * epsilon)
          } else {
            incr = incr*epsilon
          }
          lk1 = lk1+incr
        }
        current_filter[1] = current_predictive[1] * lk0
        current_filter[2] = current_predictive[2] * lk1
        l_idata = current_filter[1] + current_filter[2]
        loglikelihood_value = loglikelihood_value + log(l_idata)
        if (idata < ndata-1){
          current_filter = current_filter / l_idata
          exp_ = exp(- k * rho * gendist[idata])
          a01 = r * (1 - exp_)
          a11 = r + (1-r) * exp_
          current_predictive[2] = current_filter[1] * a01 + current_filter[2] * a11
          current_predictive[1] = 1 - current_predictive[2]
        }
    }
    return(loglikelihood_value)
}




