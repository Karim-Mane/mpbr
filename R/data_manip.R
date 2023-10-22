#' Create the `SNPdata` object
#'
#' This function generate the input data needed for whole genome SNPs data
#' genotyped from malaria parasite.
#'
#' @param vcf_file the input VCF file (required)
#' @param meta_file the metadata file (required)
#' @param output_dir the path to the folder where the output files will
#'    be stored (optional)
#' @param gaf the gene ontology annotation file (optional). If not provided, the
#'    default file obtained from the [PlasmoDB](
#' https://plasmodb.org/plasmo/app/downloads/Current_Release/Pfalciparum3D7/gaf/) # nolint: line_length_linter
#'    will be used
#' @param gff the gene annotation file (optional). If not provided, the default
#'    file obtained from the [PlasmoDB](
#' https://plasmodb.org/plasmo/app/downloads/Current_Release/Pfalciparum3D7/gff/) # nolint: line_length_linter
#'    will be used
#' @param num_threads the number of threads to be used when reading in the data
#'    from the VCF file. default is 4
#'
#' @return an object of class `SNPdata` with 5 elements
#' \enumerate{
#'   \item meta: a `data.frame` that contains the sample's metadata
#'   \item details: a `data.frame` with SNPs genomic coordinates, the fraction
#'      of missing data per SNP and the names and descriptions of the gene on
#'      which they are located.
#'   \item GT: an integer `matrix` with the genotype data. 0='reference allele',
#'      1='alternate allele', 2='mixed allele', NA='missing allele'
#'   \item vcf: the full path to the VCF file from which the data is generated
#'   \item index: an integer. this is used as some kind of label that can tell
#'      the how many times the VCF files has been modified.
#'   }
#'
#' @examples
#' \dontrun{
#'   snpdata <- get_snpdata(
#'     vcf_file   = "file.vcf.gz",
#'     meta_file  = "file.txt",
#'     output_dir = system.file("extdata", package = "mpbr")
#'  )
#' }
#' @export
#'
get_snpdata <- function(vcf_file    = NULL,
                        meta_file   = NULL,
                        output_dir  = NULL,
                        gaf         = NULL,
                        gff         = NULL,
                        num_threads = 4L) {
  checkmate::assert_file_exists(vcf_file)
  checkmate::assert_file_exists(meta_file)
  checkmate::assert_directory_exists(output_dir)
  checkmate::assert_character(gaf, any.missing = FALSE, len = 1L,
                              null.ok = TRUE)
  checkmate::assert_integer(num_threads, lower = 1L,
                            upper = (parallel::detectCores() - 1L),
                            any.missing = FALSE, null.ok = FALSE, len = 1L)
  # the user need to provide the GAF file from which the gene ontology
  # annotation will be extracted.
  # otherwise, the pre-existing GAF file will be used
  if (all(!is.null(gaf) && file.exists(gaf))) {
    go <- data.table::fread(gaf, nThread = num_threads, sep = "\t")
  } else {
    go <- data.table::fread(
      system.file("extdata", "Pf_gene_ontology.txt", package = "mpbr"),
      nThread = num_threads,
      sep     = "\t"
    )
  }

  # The user needs to provide the GFF file from which the gene names will be
  # extracted. This is converted into BED format because the GenomicRange
  # package works on such file types.
  if (all(!is.null(gff) && file.exists(gff))) {
    bed <- file.path(dirname(vcf_file), "file.bed")
    system(sprintf("gff2bed < %s > %s", gff, bed))
    bed <- data.table::fread(bed, nThread = num_threads, sep = "\t")
  } else {
    bed <- data.table::fread(
      system.file("extdata", "file.bed", package = "mpbr"),
      nThread = num_threads,
      sep = "\t"
    )
  }

  # the sample IDs will be used to create the sample metadata file.
  ids        <- file.path(output_dir, "sample_ids.txt")
  system(sprintf("bcftools query -l %s > %s", vcf_file, ids))
  sample_ids <- data.table::fread(ids, header = FALSE)

  # the genotype data will be used to create the genotype matrix and the details
  # table
  genotypes  <- file.path(output_dir, "Genotypes.txt")
  expression <- "%CHROM\t%POS\t%REF\t%ALT\t%QUAL[\t%GT]\n"
  system(sprintf("bcftools query -f'%s' %s > %s",
                 expression,
                 vcf_file,
                 genotypes))
  genotype_f <- data.table::fread(genotypes,
                                  header  = FALSE,
                                  nThread = num_threads)
  names(genotype_f) <- c("Chrom", "Pos", "Ref", "Alt", "Qual",
                         sample_ids[["V1"]])

  # the details is needed to store the genomic coordinates of the variants
  Chrom <- Pos <- Ref <- Qual <- Alt <- NULL # nolint
  details           <- genotype_f %>% dplyr::select(Chrom, Pos, Ref, Alt, Qual)
  names(sample_ids) <- "sample"
  snps              <- as.matrix(subset(genotype_f, select = -(1:5))) # nolint
  snps[snps == "0/0"]                 <- "0"
  snps[snps == "1/1"]                 <- "1"
  snps[snps == "0/1" | snps == "1/0"] <- "2"
  snps[snps == "./." | snps == ".|."] <- NA
  snps <- apply(snps, 2, function(x) as.integer(x)) # nolint
  meta <- add_metadata(sample_ids, meta_file)
  meta[["percentage_missing_sites"]]      <- colSums(is.na(snps)) / nrow(snps)
  details[["percentage_missing_samples"]] <- rowSums(is.na(snps)) / ncol(snps)

  # adding the annotation data to the details table to associate each SNPs to
  # its gene of origin together with that gene's function.
  genomic_coordinates      <- details %>% dplyr::select(Chrom, Pos)
  details[["gene"]]        <- get_gene_annotation(genomic_coordinates, go, bed)

  # we created the SNPdata class to handle easily the combined set of all the
  # data needed for downstream analyses.
  snp_table <- list(
    meta    = meta,
    details = details,
    GT      = snps,
    vcf     = vcf_file,
    index   = 0L
  )
  class(snp_table) <- "SNPdata"
  snp_table
}


#' Build the sample metadata table
#'
#' @param sample_ids a `vector` of sample IDs. This should be in the same order
#'    as they appear in the input VCF file.
#' @param metadata the path to the sample metadata file
#'
#' @return an object of class `data.frame` that contains the sample metadata
#' @keywords internal
#' @noRd
#'
add_metadata <- function(sample_ids, metadata) {
  checkmate::assert_data_frame(sample_ids, min.rows = 1L, min.cols = 1L,
                               null.ok = FALSE)
  checkmate::assert_file_exists(metadata)
  meta    <- data.table::fread(metadata, key = "sample", nThread = 4L)
  samples <- meta[["sample"]]

  # sample from the VCF file must match with those in the metadata file
  are_in_meta_file <- sample_ids[["sample"]] %in% samples
  if (!all(are_in_meta_file)) {
    warning(sprintf("Incomplete meta data - the following samples in the VCF
                    file are not found in metadata file:%s %s",
                    "\n",
                    glue::glue_collapse(
                                        sample_ids[["sample"]][!are_in_meta_file], # nolint: line_length_linter
                                        sep = ", ")),
    call. = FALSE)
  }

  # samples from the metadata file should also match with the ones from
  # the input VCF file
  are_in_vcf_file <- samples %in% sample_ids[["sample"]]
  if (!all(are_in_vcf_file)) {
    warning(sprintf("The following samples are removed from metadata file as
                    as they are not found in the VCF file: %s",
                    glue::glue_collapse(samples[!are_in_vcf_file],
                                        sep = ", ")),
            call. = FALSE)
    meta <- meta[-(!are_in_vcf_file), ]
  }

  # joining the samples IDs from the VCF file with the sample metadata
  meta <- data.frame(sample = samples) %>%
    dplyr::left_join(meta, by = "sample")

  meta
}

#' Add gene ontology and names annotation details to every SNPs in the table
#' that contains their genomic coordinates
#'
#' @param genomic_coordinates the table with SNPs genomic coordinates
#' @param go a `data.frame` with the gene ontology annotation details
#' @param bed a `data.frame` with the gene name annotation
#' @param num_cores the number of cores to be used
#'
#' @return an object of type `data.frame` with the SNPs genomic coordinates and
#'    their corresponding annotation details
#' @keywords internal
#' @noRd
#'
get_gene_annotation <- function(genomic_coordinates, go, bed, num_cores = 4L) {
  checkmate::assert_data_frame(genomic_coordinates, min.rows = 1L,
                               min.cols = 1L, null.ok = FALSE)
  checkmate::assert_data_frame(go, min.rows = 1L,
                               min.cols = 1L, null.ok = FALSE)
  checkmate::assert_data_frame(bed, min.rows = 1L,
                               min.cols = 1L, null.ok = FALSE)

  genes        <- as.character(parallel::mclapply(bed[["V10"]],
                                                  get_clean_name,
                                                  mc.cores = num_cores))
  genes        <- as.character(parallel::mclapply(genes,
                                                  rm_prf1,
                                                  mc.cores = num_cores))
  genes        <- as.character(parallel::mclapply(genes,
                                                  rm_prf2,
                                                  mc.cores = num_cores))
  genes        <- as.character(parallel::mclapply(genes,
                                                  rm_suf,
                                                  mc.cores = num_cores))
  genes        <- data.table::data.table(genes)
  genes        <- cbind(bed[["V1"]], bed[["V2"]], bed[["V3"]], genes)
  names(genes) <- c("chrom", "start", "end", "gene_id")
  go           <- subset(go, select = c(2, 10)) # nolint
  names(go)    <- c("gene_id", "gene_name")
  data.table::setkey(go, "gene_id")
  data.table::setkey(genes, "gene_id")

  chrom <- gene_id <- gene_name <- NULL
  test         <- genes %>% dplyr::left_join(go)
  test         <- dplyr::distinct(test, chrom, start, end, gene_id, gene_name)
  resultat     <- gene_annotation(test, genomic_coordinates)
  resultat     <- gsub("NA:", "", resultat, fixed = TRUE)
  resultat
}


#' clean gene names
#'
#' @param y the gene name
#'
#' @return a string with the cleaned gene name
#' @keywords internal
#' @noRd
#'
get_clean_name <- function(y) {
  checkmate::assert_character(y, any.missing = FALSE, null.ok = FALSE)
  unlist(strsplit(
                  unlist(strsplit(y, "ID=", fixed = TRUE))[[2L]],
                  ";", fixed = TRUE))[[1L]]
}

#' remove the "exon_" prefix from character string
#'
#' @param x the input character string
#'
#' @return the input character string without the prefix "exon_"
#' @keywords internal
#' @noRd
#'
rm_prf1 <- function(x) {
  checkmate::assert_character(x, any.missing = FALSE, null.ok = FALSE)
  as.character(gsub("exon_", "", x, fixed = TRUE))
}

#' remove the "utr_" preffix from a character
#'
#' @param x the input character string
#'
#' @return the input character string without the preffix "utr_"
#' @keywords internal
#' @noRd
#'
rm_prf2 <- function(x) {
  checkmate::assert_character(x, any.missing = FALSE, null.ok = FALSE)
  as.character(gsub("utr_", "", x, fixed = TRUE))
}

#' split a character on '.'
#'
#' @param x the input character string
#'
#' @return a `vector` of `character` obtained after splitting the input
#'    character based on '.'
#' @keywords internal
#' @noRd
#'
rm_suf <- function(x) {
  checkmate::assert_character(x, any.missing = FALSE, null.ok = FALSE)
  as.character(unlist(strsplit(x, ".", fixed = TRUE))[[1L]])
}

#' Associate each SNPs to a gene name on which it belongs
#'
#' @param target_gtf a `data.frame` with the gene annotation
#' @param genomic_coordinates a `data.frame` with the SNPs genomic coordinates
#'
#' @return a `vector` with gene annotation. This should be of the same length as
#'    the number of SNPs.
#' @keywords internal
#' @noRd
#'
gene_annotation <- function(target_gtf, genomic_coordinates) {
  checkmate::assert_data_frame(target_gtf, min.rows = 1L, min.cols = 1L,
                               null.ok = FALSE)
  checkmate::assert_data_frame(genomic_coordinates,
                               min.rows = 1L,
                               min.cols = 1L,
                               null.ok  = FALSE)

  names(genomic_coordinates)   <- c("chrom", "start")
  genomic_coordinates[["end"]] <- genomic_coordinates[["start"]]
  subject                      <- IRanges::IRanges(target_gtf[["start"]],
                                                   target_gtf[["end"]])
  query       <- IRanges::IRanges(genomic_coordinates[["start"]],
                                  genomic_coordinates[["end"]])
  my_overlaps <-
    data.table::data.table(as.matrix(GenomicRanges::findOverlaps(
                                                                 query,
                                                                 subject,
                                                                 type = "within"))) #nolint: line_length_linter
  my_overlaps[["gene"]] <- target_gtf[["gene_name"]][my_overlaps[["subjectHits"]]] #nolint: line_length_linter
  gene <- queryHits <- NULL # nolint
  the_genes <- my_overlaps[, paste(unique(gene), collapse = ":"),
                           by = queryHits]
  names(the_genes)[[2L]]           <- "gene"
  genomic_coordinates[["gene"]]    <- NA
  genomic_coordinates[["gene"]][the_genes[["queryHits"]]] <- the_genes[["gene"]]

  genomic_coordinates[["gene"]]
}

#' Print the `SNPdata` object
#'
#' @param snpdata the `SNPdata` object
#'
#' @returns prints the structure of the `SNPdata` object
#'
#' @examples
#' \dontrun{
#'   print(snpdata)
#' }
#'
#' @export
print.SNPdata <- function(snpdata) {
  checkmate::assert_class(snpdata, "SNPdata", null.ok = FALSE)
  print(head(snpdata[["meta"]]))
  print(head(snpdata[["details"]]))
  message(sprintf("\nData contains: %d samples for %d snp loci\n",
                  dim(snpdata[["GT"]])[[2L]],
                  dim(snpdata[["GT"]])[[1L]]))
  message(sprintf("\nData is generated from: %s\n", snpdata[["vcf"]]))
}

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
#'   maf_cutoff          = 0.01
#'  )
#'  }
#' @export
#'
filter_snps_samples <- function(snpdata, min_qual   = 10L,
                                max_missing_sites   = 0.2,
                                max_missing_samples = 0.2,
                                maf_cutoff          = 0.01) {
  checkmate::assert_class(snpdata, "SNPdata", null.ok = FALSE)
  checkmate::assert_integer(min_qual, lower = 10L, any.missing = FALSE,
                            null.ok = FALSE)
  checkmate::assert_numeric(max_missing_sites, lower = 0L, upper = 1L,
                            finite = TRUE, any.missing = FALSE, null.ok = FALSE,
                            len = 1L)
  checkmate::assert_numeric(max_missing_samples, lower = 0L, upper = 1L,
                            finite = TRUE, any.missing = FALSE, null.ok = FALSE,
                            len = 1L)
  checkmate::assert_numeric(maf_cutoff, lower = 0L, upper = 1L,
                            finite = TRUE, any.missing = FALSE, null.ok = FALSE,
                            len = 1L)
  x      <- snpdata[["details"]]
  fields <- c("GT", "Phased", "Phased_Imputed")
  if (all(missing(min_qual) &&
            missing(max_missing_sites) &&
            missing(max_missing_samples))) {
    return(snpdata)
  } else {
    idx <- which(x[["Qual"]] >= min_qual &
                   x[["percentage_missing_samples"]] <= max_missing_samples &
                   x[["MAF"]] >= maf_cutoff)
    if (all(length(idx) > 0L && length(idx) < nrow(snpdata[["details"]]))) {
      x <- x[idx, ]
      snpdata[["details"]] <- x
      for (field in fields) {
        if (field %in% names(snpdata)) {
          snpdata[[field]] <- snpdata[[field]][idx, ]
        }
      }
      Chrom <- Pos <- NULL # nolint
      f2c        <- x %>% dplyr::select(Chrom, Pos)
      output_dir <- dirname(snpdata[["vcf"]])
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
    } else if (length(idx) == 0L) {
      stop("\nNo locus in VCF file has satisfied specified the QC metrics")
    } else if (length(idx) == nrow(snpdata[["details"]])) {
      message("all loci have satisfied the specified QC metrics")
    }

    idx <- which(snpdata[["meta"]][["percentage_missing_sites"]] <=
                   max_missing_sites)
    if (all(length(idx) > 0L & length(idx) < nrow(snpdata[["meta"]]))) {
      message("\nthe following samples will be removed:\n",
              paste(snpdata[["meta"]][["sample"]], collapse = "\n"))
      snpdata[["meta"]] <- snpdata[["meta"]][idx, ]
      data.table::fwrite(snpdata[["meta"]][["sample"]],
                         file.path(output_dir, "samples_to_be_dropped.txt"),
                         col.names = FALSE, row.names = FALSE,
                         quote = FALSE, sep = "\t", nThread = 4L)
      snpdata[["vcf"]] <- remove_samples_from_vcf(snpdata[["vcf"]],
                                                  "samples_to_be_dropped.txt",
                                                  output_dir,
                                                  index = snpdata[["index"]])
    } else if (length(idx) == 0L) {
      stop("\nNo sample in VCF file has satisfied the specified QC metrics.")
    } else if (length(idx) == nrow(snpdata[["meta"]])) {
      message("\nAll samples have satisfied the specified QC metrics.")
    }
  }

  snpdata[["index"]] <- snpdata[["index"]] + 1L
  snpdata
}

#' Calculate minor allele frequency (MAF)
#'
#' Uses the `SNPdata` object to calculate the MAF at every loci
#'
#' @param snpdata a `SNPdata` object
#' @param include_het whether to account for the heterozygous allele or not.
#'    This is only used when `mat_name = "GT"`
#' @param mat_name the name of the matrix to use. default is "GT"
#'
#' @return a `SNPdata` object with 2 additional columns in the **details**
#'    table.
#' \enumerate{
#'   \item MAF: minor allele frequency of each SNPs
#'   \item MAF_allele: 1 if the alternate allele is the minor allele.
#'         0 otherwise
#' }
#'
#' @details if `include_het = FALSE`, the mixed allele will not be considered in
#'    the MAF calculation
#'
#' @examples
#' \dontrun{
#'   snpdata <- compute_maf(
#'    snpdata,
#'    include_het = FALSE,
#'    mat_name    = "GT"
#'  )
#'  }
#' @export
#'
compute_maf <- function(snpdata, include_het = FALSE, mat_name = "GT") {
  checkmate::assert_class(snpdata, "SNPdata", null.ok = FALSE)
  checkmate::assert_logical(include_het, any.missing = FALSE, len = 1L,
                            null.ok = FALSE)
  checkmate::assert_character(mat_name, any.missing = FALSE, null.ok = FALSE,
                              len = 1L)
  x   <- snpdata[[mat_name]]
  ref <- rowSums(x == 0L, na.rm = TRUE)
  alt <- rowSums(x == 1L, na.rm = TRUE)
  het <- rowSums(x == 2L, na.rm = TRUE)
  tmp_mat <- ifelse(!include_het, cbind(ref, alt), cbind(ref, alt, het))
  res <- apply(tmp_mat, 1L, get_maf)
  if (!("MAF" %in% names(snpdata[["details"]]))) {
    snpdata[["details"]][["MAF"]]           <- as.numeric(res[1L, ])
    snpdata[["details"]][["MAF_allele"]] <-
      as.factor(as.character(as.numeric(round(res[2L, ]))))
    levels(snpdata[["details"]][["MAF_allele"]]) <-
      dplyr::recode_factor(snpdata[["details"]][["MAF_allele"]], REF = "0",
                           ALT = "1", HET = "2", REF_ALT = "3",
                           REF_ALT_HET = "4")
  } else {
    new_maf <- paste0("MAF_", mat_name)
    snpdata[["details"]][[new_maf]] <- as.numeric(res[1L, ])
  }

  snpdata
}

#' get the minor allele frequency (MAF) and the corresponding allele
#'
#' @param mat a matrix with the 2 or 3 columns. Every column should contain the
#'    the count of the allele of interest across all samples.
#'
#' @return a `vector` of the MAF and the corresponding allele
#' @keywords internal
#' @noRd
#'
get_maf <- function(mat) {
  checkmate::assert_matrix(mat, mode = "numeric", min.rows = 1L, min.cols = 2L,
                           null.ok = FALSE)
  if (length(mat) == 2L) {
    if (mat[[1L]] < mat[[2L]]) {
      maf    <- mat[[1L]] / sum(mat[[1L]], mat[[2L]])
      allele <- 0L
    } else if (mat[[1L]] > mat[[2L]]) {
      maf    <- mat[[2L]] / sum(mat[[1L]], mat[[2L]])
      allele <- 1L
    } else {
      maf    <- mat[[2L]] / sum(mat[[1L]], mat[[2L]])
      allele <- 3L
    }
  } else {
    if (mat[[1L]] < mat[[2L]]) {
      minor  <- mat[[1L]]
      allele <- 0L
    } else if (mat[[1L]] > mat[[2L]]) {
      minor  <- mat[[2L]]
      allele <- 1L
    } else {
      minor  <- mat[[2L]]
      allele <- 3L
    }

    if (minor < mat[[3L]]) {
      maf    <- minor / sum(mat[[1L]], mat[[2L]], mat[[3L]])
      allele <- 3L
    } else if (minor > mat[[3L]]) {
      maf    <- mat[[3L]] / sum(mat[[1L]], mat[[2L]], mat[[3L]])
      allele <- 2L
    } else {
      maf    <- mat[[3L]] / sum(mat[[1L]], mat[[2L]], mat[[3L]])
      allele <- 4L
    }
  }

  return(c(maf, allele))
}


#' Calculate the complexity of the infection (COI) in every sample
#'
#' The COI is estimated here through the within host genetic diversity (Fws).
#' using the {moimix} R package.
#'
#' @param snpdata a `SNPdata` object
#'
#' @return a `SNPdata` object with 2 additional columns in the meta table
#' \enumerate{
#'   \item Fws: within host genetic diversity value
#'   \item COI: the complexity of infection: 1 for Fws>0.95, 2 for Fws<=0.95
#' }
#'
#' @examples
#' \dontrun{
#'   snpdata <- calculate_fws(snpdata)
#'  }
#'
#' @export
#'
calculate_fws <- function(snpdata) {
  checkmate::assert_class(snpdata, "SNPdata", null.ok = FALSE)
  vcf        <- snpdata[["vcf"]]
  gds_file   <- file.path(dirname(vcf), "data.gds")
  SeqArray::seqVCF2GDS(vcf, gds_file)
  my_vcf    <- SeqArray::seqOpen(gds_file)

  # filter out the SNPs from the apicoplast and the mitochondrial chromosomes
  sample_id <- SeqArray::seqGetData(my_vcf, "sample_id")
  coords    <- moimix::getCoordinates(my_vcf)
  SeqArray::seqSetFilter(my_vcf,
                         variant.id = coords[["variant.id"]][coords[["chromosome"]] != "Pf3D7_API_v3"]) # nolint: line_length_linter
  SeqArray::seqSetFilter(my_vcf,
                         variant.id = coords[["variant.id"]][coords[["chromosome"]] != "Pf_M76611"]) # nolint: line_length_linter

  # calculate the Fws and the COI
  fws_overall        <- moimix::getFws(my_vcf)
  fws_overall        <- data.table::data.table(cbind(as.character(sample_id),
                                                     as.numeric(fws_overall)))
  names(fws_overall) <- c("sample", "Fws")
  meta <- data.frame(snpdata[["meta"]]) %>%
    dplyr::left_join(fws_overall, by = "sample")
  meta[["COI"]]      <- 1L
  meta[which(meta[["Fws"]] <= 0.95), ][["COI"]] <- 2L
  snpdata[["meta"]]  <- meta
  snpdata
}

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
  checkmate::assert_integer(nsim, lower = 1L, any.missing = FALSE,
                            null.ok = FALSE, len = 1L)
  vcf          <- snpdata[["vcf"]]
  expression   <- '%CHROM\t%POS[\t%AD]\n' # nolint
  tmp          <- file.path(dirname(vcf), "tmp")
  system(sprintf("mkdir -p %s", tmp))
  ad           <- file.path(tmp, "AllelicDepth.txt")
  system(sprintf("bcftools query -f'%s' %s > %s", expression, vcf, ad))
  depth        <- data.table::fread(ad, nThread = 4L)
  depth        <- as.matrix(subset(depth, select = -c(1L:2L))) # nolint
  path         <- file.path(dirname(vcf), "phasing")
  system(sprintf("mkdir -p %s", path))
  correlations <- numeric(length = nsim)
  pb           <- txtProgressBar(min = 0L, max = nsim, initial = 0L,
                                 style = 3L, char = "*")
  for (i in 1L:nsim) {
    tmp_snpdata <- snpdata
    mat         <- apply(tmp_snpdata[["GT"]], 1L, phase_data, depth = depth)
    tmp_snpdata[["Phased"]] <- t(mat)
    saveRDS(t(mat), file.path(path, paste0("sim", i, ".RDS")))
    res_snpdata <- compute_maf(tmp_snpdata,
                               include_het = FALSE,
                               mat_name = "Phased")
    correlations[i] <- cor(res_snpdata[["details"]][["MAF_Phased"]],
                           res_snpdata[["details"]][["MAF"]])
    setTxtProgressBar(pb, i)
  }
  close(pb)
  idx <- which(correlations == max(correlations, na.rm = TRUE))
  snpdata[["Phased"]] <- readRDS(file.path(path,
                                           paste0("sim", idx[[1L]], ".RDS")))
  system(sprintf("rm -rf %s", path))
  snpdata
}

#' Phase the mixed genotypes
#'
#' @param genotype a vector of genotype data
#' @param depth a vector of allelic depth
#'
#' @return a vector of phased genotypes
#' @keywords internal
#' @noRd
#'
phase_data <- function(genotype, depth) {
  checkmate::assert_vector(genotype, min.len = 1L, null.ok = FALSE)
  checkmate::assert_vector(depth, min.len = 1L, null.ok = FALSE)
  idx <- as.numeric(which(genotype == 2L))

  for (j in idx) {
    ref <- as.numeric(unlist(strsplit(depth[j], ',', fixed = TRUE))[[1L]]) # nolint: quotes_linters
    alt <- as.numeric(unlist(strsplit(depth[j], ',', fixed = TRUE))[[2L]]) # nolint: quotes_linters
    if (all(ref == 0L && alt == 0L)) {
      ref_count <- sum(genotype == 0L, na.rm = TRUE)
      alt_count <- sum(genotype == 1L, na.rm = TRUE)
      if (ref_count < alt_count) {
        genotype[j] <- 0L
      } else if (ref_count > alt_count) {
        genotype[j] <- 1L
      } else {
        genotype[j] <- statip::rbern(1L, ref_count / (ref_count + alt_count))
      }
    } else if (all(ref != 0L && alt != 0L)) {
      if ((ref + alt) >= 5L && (ref >= (2L * alt) || alt >= (2L * ref))) {
        if (ref < alt) {
          genotype[j] <- 0L
        } else if (ref > alt) {
          genotype[j] <- 1L
        } else {
          ref_count <- sum(genotype == 0L, na.rm = TRUE)
          alt_count <- sum(genotype == 1L, na.rm = TRUE)
          if (ref_count < alt_count) {
            genotype[j] <- 0L
          } else if (ref_count > alt_count) {
            genotype[j] <- 1L
          } else {
            genotype[j] <- statip::rbern(1L,
                                         ref_count / (ref_count + alt_count))
          }
        }
      } else {
        ref_count <- sum(genotype == 0L, na.rm = TRUE)
        alt_count <- sum(genotype == 1L, na.rm = TRUE)
        if (ref_count < alt_count) {
          genotype[j] <- 0L
        } else if (ref_count > alt_count) {
          genotype[j] <- 1L
        } else {
          genotype[j] <- statip::rbern(1L,
                                       ref_count / (ref_count + alt_count))
        }
      }
    } else if (all(ref == 0L || alt == 0L)) {
      ref_count <- sum(genotype == 0L, na.rm = TRUE)
      alt_count <- sum(genotype == 1L, na.rm = TRUE)
      if (all(ref == 0L && alt >= 5L)) {
        genotype[j] <- 1L
      } else if (all(ref == 0L && alt < 5L)) {
        genotype[j] <- statip::rbern(1L,
                                     alt_count / (ref_count + alt_count))
      }
      if (all(alt == 0L && ref >= 5L)) {
        genotype[j] <- 0L
      } else if (all(alt == 0L && ref < 5L)) {
        genotype[j] <- statip::rbern(1L, ref_count / (ref_count + alt_count))
      }
    }
  }
  genotype
}


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
  checkmate::assert_integer(nsim, lower = 1L, any.missing = FALSE,
                            null.ok = FALSE, len = 1L)
  message("The missing genotypes will be imputed from ", genotype, " table.\n")
  field <- genotype
  path  <- file.path(dirname(snpdata[["vcf"]]), "imputing")
  system(sprintf("mkdir -p %s", path))
  correlations <- numeric(length = nsim)
  pb <- txtProgressBar(min = 0L, max = nsim, initial = 0L,
                       style = 3L, char = "*")
  for (i in 1L:nsim) {
    tmp_snpdata <- snpdata
    mat         <- apply(tmp_snpdata[[field]], 1L, impute)
    tmp_snpdata[["Imputed"]] <- t(mat)
    saveRDS(t(mat), file.path(path, paste0("sim", i, ".RDS")))
    res_snpdata <- compute_maf(tmp_snpdata, include_het = FALSE,
                               mat_name = "Imputed")
    correlations[i] <- cor(res_snpdata[["details"]][["MAF_Imputed"]],
                           res_snpdata[["details"]][["MAF"]])
    setTxtProgressBar(pb, i)
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
#' @noRd
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

#' Select data from specified chromosomes
#'
#' @param snpdata a `SNPdata` object
#' @param chrom a vector of chromosome names
#'
#' @return a `SNPdata` object with only the data from the specified chromosomes
#'
#' @examples
#' \dontrun{
#'   chrom_snpdata <- select_chrom(snpdata, chrom = "Pf3D7_07_v3")
#'  }
#' @export
#'
select_chrom <- function(snpdata, chrom = "all") {
  checkmate::assert_class(snpdata, "SNPdata", null.ok = FALSE)
  checkmate::assert_vector(chrom, any.missing = FALSE, null.ok = FALSE,
                           min.len = 1L)
  system(sprintf("tabix %s", snpdata[["vcf"]]))
  m      <- which(names(snpdata) %in% c("meta", "vcf", "index"))
  fields <- names(snpdata)[-m]
  if (chrom == "all") {
    return(snpdata)
  }
  res    <- list()
  for (chr in chrom) {
    chrom_snpdata <- snpdata
    idx           <- which(chrom_snpdata[["details"]][["Chrom"]] == chr)
    for (field in fields) {
      res[[field]] <- rbind(res[[field]], chrom_snpdata[[field]][idx, ])
    }
  }
  chrom_vcf <- file.path(dirname(chrom_snpdata[["vcf"]]), "target_chrom_vcf.gz")
  if (file.exists(chrom_vcf)) {
    system(sprintf("rm -f %s", chrom_vcf))
  }
  if (length(chrom) > 1L) {
    tmp_xme <- file.path(dirname(chrom_snpdata[["vcf"]]), "target_chrom.txt")
    data.table::fwrite(chrom, tmp_xme, col.names = FALSE,
                       row.names = FALSE, quote = FALSE, sep = "\t")
    system(sprintf("bcftools view -R %s %s -o %s -O z",
                   tmp_xme, chrom_snpdata[["vcf"]], chrom_vcf))
  } else {
    system(sprintf("bcftools view -r\"%s\" %s -o %s -O z",
                   chrom, chrom_snpdata[["vcf"]], chrom_vcf))
  }
  res[["vcf"]]        <- chrom_vcf
  snpdata[["meta"]]   <- snpdata[["meta"]]
  class(res)          <- "SNPdata"
  res
}

#' Drop a of set SNPs from a `SNPdata` object
#'
#' @param snpdata a `SNPdata` object
#' @param snp_to_be_dropped a `data.frame` with 2 columns" "Chrom" and "Pos"
#' @param chrom the chromosome from which loci should be dropped
#' @param start the starting position of the region to be discarded
#' @param end the end position of the region to be discarded
#'
#' @return a `SNPdata` object where the specified SNPs have been removed
#'
#' @examples
#' \dontrun{
#'   snpdata <- drop_snps(
#'     snpdata,
#'     snp_to_be_dropped = NULL,
#'     chrom             = "Pf3D7_05_v3",
#'     start             = 100,
#'     end               = 500
#'   )
#'  }
#' @details when snp_to_be_dropped is not set to NA (i.e. the genomic
#'     coordinates of snps to be removed are in a data frame), then the rest of
#'     the arguments can be ignored or set to NA (chrom=NA, start=NA, end=NA)
#' @export
#'
drop_snps <- function(snpdata, snp_to_be_dropped = NULL,
                      chrom = NULL, start = NULL, end = NULL) {
  checkmate::assert_class(snpdata, "SNPdata", null.ok = FALSE)
  checkmate::assert_data_frame(snp_to_be_dropped, ncols = 2L, null.ok = TRUE,
                               col.names = c("Chrom", "Pos"))
  checkmate::assert_vector(chrom, min.len = 1L, any.missing = FALSE,
                           null.ok = TRUE)
  checkmate::assert_vector(start, min.len = 1L, any.missing = FALSE,
                           null.ok = TRUE)
  checkmate::assert_vector(end, min.len = 1L, any.missing = FALSE,
                           null.ok = TRUE)
  if (all(is.null(snp_to_be_dropped) &&
            is.null(chrom) &&
            is.null(start) &&
            is.null(end))) {
    stop("Please provide genomic coordinates of loci to be removed.")
  }
  if (!is.null(snp_to_be_dropped) &&
        (all(is.null(chrom) && is.null(start) && is.null(end)))) {
    idx <- which(snpdata[["details"]][["Chrom"]] %in%
                   snp_to_be_dropped[["Chrom"]] &
                   snpdata[["details"]][["Pos"]] %in%
                     snp_to_be_dropped[["Pos"]])
    m      <- which(names(snpdata) %in% c("meta", "vcf", "index"))
    fields <- names(snpdata)[-m]
    for (field in fields) {
      tmp              <- snpdata[[field]][-idx, ]
      snpdata[[field]] <- tmp
    }
    Chrom <- Pos <- NULL # nolint: object_name_linter
    f2c      <- snpdata[["details"]] %>%
      dplyr::select(Chrom, Pos)
    tmp_file <- file.path(dirname(snpdata[["vcf"]]), "tmp.txt")
    data.table::fwrite(f2c, tmp_file, col.names = FALSE, row.names = FALSE,
                       quote = FALSE, sep = "\t", nThread = 4L)
    snpdata[["vcf"]] <- remove_snps_from_vcf(snpdata[["vcf"]], "tmp.txt",
                                             path = dirname(snpdata[["vcf"]]),
                                             index = snpdata[["index"]])
  } else if (all(!is.null(chrom) && !is.null(start) &&
                   !is.null(end)) &&
               is.na(snp_to_be_dropped)) {
    idx <- which(snpdata[["details"]][["Chrom"]] == chrom &
                   (snpdata[["details"]][["Pos"]] >= start &
                      snpdata[["details"]][["Pos"]] <= end))
    if (length(idx) > 0L) {
      m      <- which(names(snpdata) %in% c("meta", "vcf", "index"))
      fields <- names(snpdata)[-m]
      for (field in fields) {
        tmp              <- snpdata[[field]][-idx, ]
        snpdata[[field]] <- tmp
      }
      f2c      <- snpdata[["details"]] %>%
        dplyr::select(Chrom, Pos)
      tmp_file <- file.path(dirname(snpdata[["vcf"]]), "tmp.txt")
      data.table::fwrite(f2c, tmp_file, col.names = FALSE,
                         row.names = FALSE, quote = FALSE,
                         sep = "\t", nThread = 4L)
      snpdata[["vcf"]] <- remove_snps_from_vcf(snpdata[["vcf"]],
                                               "tmp.txt",
                                               path  = dirname(snpdata[["vcf"]]), # nolint: line_length_linter
                                               index = snpdata[["index"]])
      system(sprintf("rm -f %s", tmp_file))
      message("\n", length(idx), "loci have been successfully removed")
    } else {
      stop("There is no loci overlapping the specified region")
    }
  } else {
    stop("Incorrect genomics coordinates")
  }
  snpdata[["index"]] <- snpdata[["index"]] + 1L
  snpdata
}

#' Delete a set of SNPs from a VCF file
#'
#' @param vcf the input VCF file
#' @param loci_to_be_retained the name of the file that contains the set of loci
#'    to be retained.
#' @param path the path to the directory that contains the input files. This is
#'    also where the output files will be stored.
#' @param index the index of the input VCF file
#'
#' @return the path to the filtered VCF file
#' @keywords internal
#' @noRd
#'
remove_snps_from_vcf <- function(vcf, loci_to_be_retained, path, index = 1L) {
  checkmate::assert_file_exists(vcf, extension = c(".vcf", "vcf.gz"))
  checkmate::assert_character(loci_to_be_retained, any.missing = FALSE,
                              null.ok = FALSE)
  checkmate::assert_directory_exists(path)
  checkmate::assert_integer(index, lower = 0L, len = 1L, null.ok = FALSE,
                            any.missing = FALSE)
  target_loci <- file.path(path, loci_to_be_retained)
  header      <- file.path(path, "Header.txt")
  body        <- file.path(path, "Body.txt")
  correct_rows <- file.path(path, "Good_snps.txt")
  filtered_vcf <- file.path(path, paste0("Filtered_snps_", index, ".vcf"))
  system(sprintf("bcftools view -h %s > %s", vcf, header))
  system(sprintf("bcftools view -H %s > %s", vcf, body))
  system(sprintf("awk -F'\t' 'NR==FNR{c[$1$2]++;next};c[$1$2] > 0' %s %s > %s",
                 target_loci, body, correct_rows))
  system(sprintf("cat %s %s > %s", header, correct_rows, filtered_vcf))
  system(sprintf("bgzip %s", filtered_vcf))
  system(sprintf("rm -f %s %s %s", header, body, correct_rows))
  filtered_vcf <- file.path(path, paste0("Filtered_snps_", index, ".vcf.gz"))

  return(as.character(filtered_vcf))
}

#' Drop a specific set of samples from the `SNPdata` object
#'
#' @param snpdata a `SNPdata` object
#' @param samples_to_be_dropped a vector of samples to be dropped
#'
#' @return a `SNPdata` object where the specified samples have been removed
#' @examples
#' \dontrun{
#'   snpdata <- drop_samples(snpdata, samples_to_be_dropped)
#'  }
#'
#' @export
drop_samples <- function(snpdata, samples_to_be_dropped) {
  checkmate::assert_class(snpdata, "SNPdata", null.ok = FALSE)
  checkmate::assert_vector(samples_to_be_dropped, any.missing = FALSE,
                           null.ok = FALSE, min.len = 1L)
  if (!all(samples_to_be_dropped %in% snpdata[["meta"]][["sample"]])) {
    stop("Some samples in the provided vector are not found in the current
         data")
  }
  idx             <- match(samples_to_be_dropped, snpdata[["meta"]][["sample"]])
  tmp_meta          <- snpdata[["meta"]]
  tmp_meta          <- tmp_meta[-(idx), ]
  snpdata[["meta"]] <- tmp_meta
  fields            <- c("GT", "Phased", "Imputed")
  for (field in fields) {
    if (field %in% names(snpdata)) {
      idx              <- match(samples_to_be_dropped,
                                colnames(snpdata[[field]]))
      m                <- seq_len(ncol(snpdata[[field]]))
      m                <- m[-idx]
      tmp_meta         <- snpdata[[field]]
      tmp_meta         <- tmp_meta[, m]
      snpdata[[field]] <- tmp_meta
    }
  }
  tmp_file <- file.path(dirname(snpdata[["vcf"]]), "tmp.txt")
  write.table(snpdata[["meta"]][["sample"]], tmp_file, col.names = FALSE,
              row.names = FALSE, quote = FALSE, sep = "\t")
  snpdata[["vcf"]]   <- remove_samples_from_vcf(snpdata[["vcf"]], "tmp.txt",
                                                path  = dirname(snpdata[["vcf"]]), # nolint: line_length_linter
                                                index = snpdata[["index"]])
  snpdata[["index"]] <- snpdata[["index"]] + 1L
  snpdata
}

#' Remove samples from VCF file
#'
#' @param vcf the input VCF file
#' @param samples_to_be_retained a vector of samples to be retained
#' @param path the path folder that contains the input data. This will be also
#'    the folder where the output files will be stored.
#' @param index the index of the input VCF file.
#'
#' @return the path to filtered VCF file
#' @keywords internal
#' @noRd
#'
remove_samples_from_vcf <- function(vcf, samples_to_be_retained,
                                    path, index = 1L) {
  checkmate::assert_file_exists(vcf)
  checkmate::assert_vector(samples_to_be_retained, any.missing = FALSE,
                           null.ok = FALSE, min.len = 1L)
  checkmate::assert_directory_exists(path)
  checkmate::assert_integer(index, lower = 0L, any.missing = FALSE,
                            null.ok = FALSE, len = 1L)
  target_samples <- file.path(path, samples_to_be_retained)
  post_qc        <- file.path(path, paste0("Post_QC_", index, ".vcf.gz"))
  system(sprintf("bcftools view -S %s %s -o %s -O z",
                 target_samples, vcf, post_qc))
  system(sprintf("tabix %s", post_qc))
  as.character(post_qc)
}
