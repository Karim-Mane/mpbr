#' Create the `SNPdata` object
#'
#' This function generate the input data needed for whole genome SNPs data
#' genotyped from malaria parasite.
#'
#' @param vcf_file the input VCF file (required)
#' @param meta_file the metadata file (required)
#' @param output_dir the path to the folder where the output files will
#'    be stored (optional)
#' @param gof the gene ontology annotation file (optional). If not provided, the
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
#'     vcf_file   = system.file("extdata", "Input_Data.vcf.gz",
#'                              package = "mpbr"),
#'     meta_file  = system.file("extdata", "SampleMetadata.RDS",
#'                              package = "mpbr"),
#'     output_dir = tempdir()
#'  )
#' }
#' @export
#' @importFrom magrittr %>%
#'
get_snpdata <- function(vcf_file    = NULL,
                        meta_file   = NULL,
                        output_dir  = NULL,
                        gof         = NULL,
                        gff         = NULL,
                        num_threads = 4L) {
  checkmate::assert_file_exists(vcf_file)
  checkmate::assert_file_exists(meta_file)
  checkmate::assert_directory_exists(output_dir)
  checkmate::assert_character(gof, any.missing = FALSE, len = 1L,
                              null.ok = TRUE)
  checkmate::assert_numeric(num_threads, lower = 1L,
                            upper = (parallel::detectCores() - 1L),
                            any.missing = FALSE, null.ok = FALSE, len = 1L)
  # the user need to provide the GOF file from which the gene ontology
  # annotation will be extracted.
  # otherwise, the pre-existing GOF file will be used
  if (!is.null(gof) && file.exists(gof)) {
    if (grepl(".RDS", basename(gof))) {
      go <- readRDS(gof)
    } else {
      data.table::fread(gof, nThread = num_threads, sep = "\t")
    }
  } else {
    go <- readRDS(system.file("extdata", "pf_gene_ontology.RDS",
                              package = "mpbr"))
  }
  
  # The user needs to provide the GFF file from which the gene names will be
  # extracted. This is converted into BED format because the `GenomicRange`
  # package requires that file type.
  if (!is.null(gff) && file.exists(gff)) {
    bed <- file.path(output_dir, "file.bed")
    system(sprintf("gff2bed < %s > %s", gff, bed))
    bed <- data.table::fread(bed, nThread = num_threads, sep = "\t")
  } else {
    bed <- readRDS(system.file("extdata", "PlasmoDB-56_Pfalciparum3D7.RDS",
                               package = "mpbr"))
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
  Chrom <- Pos <- Ref <- Qual <- Alt <- NULL # nolint: object_name_linter
  details           <- genotype_f %>% dplyr::select(Chrom, Pos, Ref, Alt, Qual)
  names(sample_ids) <- "sample"
  snps              <- as.matrix(subset(genotype_f, select = -(1L:5L)))
  snps[snps == "0/0"]                 <- "0"
  snps[snps == "1/1"]                 <- "1"
  snps[snps == "0/1" | snps == "1/0"] <- "2"
  snps[snps == "./." | snps == ".|."] <- NA
  snps <- apply(snps, 2L, as.integer)
  meta <- add_metadata(sample_ids, meta_file, output_dir)
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