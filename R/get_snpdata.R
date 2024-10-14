#' Create the `SNPdata` object
#'
#' The function generates the input data needed for population genetics
#' analyses using whole genome SNPs data genotyped from the malaria parasite.
#'
#' @param vcf_file The path to the input VCF file (required)
#' @param meta_file The path to the sample metadata file (required)
#' @param output_dir The path to the folder where the output and temporary files
#'    will be stored (optional)
#' @param gof The gene ontology annotation file (optional). If not provided, the
#'    default file obtained from the [PlasmoDB](
#' https://plasmodb.org/plasmo/app/downloads/Current_Release/Pfalciparum3D7/gaf/) # nolint: line_length_linter
#'    will be used
#' @param gff The gene annotation file (optional). If not provided, the default
#'    file obtained from the [PlasmoDB](
#' https://plasmodb.org/plasmo/app/downloads/Current_Release/Pfalciparum3D7/gff/) # nolint: line_length_linter
#'    will be used
#' @param num_threads The number of threads to be used when reading in the data
#'    from the VCF file. default is 4
#'
#' @return An object of class `SNPdata` with the following 4 elements:
#' \enumerate{
#'   \item meta: A `data.frame` that contains the sample's metadata
#'   \item details: A `data.frame` with the SNPs genomic coordinates, the
#'      fraction of missing data per SNP and the names and descriptions of the
#'      gene on which they are located.
#'   \item GT: An integer `matrix` with the genotype data. 0='reference allele',
#'      1='alternate allele', 2='mixed allele', NA='missing allele'
#'   \item vcf: the full path to the VCF file from which the data is generated.
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
    if (grepl(".RDS", basename(gff))) {
      # read annotation from existing file
      bed <- readRDS(gff)
    }
    if (grepl(".gz", basename(gff))) {
      # create bed file from downloaded annotation file
      bed <- file.path(output_dir, "file.bed")
      system(sprintf("gff2bed < %s > %s", gff, bed))
      bed <- data.table::fread(bed, nThread = num_threads, sep = "\t")
    }
  } else {
    # read annotation from existing file
    bed <- readRDS(
      system.file(
        "extdata",
        "PlasmoDB-56_Pfalciparum3D7.RDS",
        package = "mpbr"
      )
    )
  }

  # the sample IDs will be used to create the sample metadata file.
  # the linux commands are different from the windows commands
  # some part of the script will be tailored to the user's operating system
  os_type <- Sys.info()["sysname"]

  if (os_type == "Windows") {
    # use gzip -dc for Windows
    command <- sprintf("gzip -dc %s | grep '^#CHROM' | cut -f 10-", vcf_file)
  } else {
    # use gzcat for Unix-like systems (macOS, Linux)
    command <- sprintf(
      "gzcat %s | grep '^#CHROM' | cut -d$'\t' -f 10-",
      vcf_file
    )
  }
  sample_ids <- as.character(data.table::fread(
    cmd = command,
    header = FALSE
  ))

  # the genotype data will be used to create the genotype matrix and the details
  # table
  genotype_data <- extract_genotype(vcf_file)
  names(genotype_data) <- c("Chrom", "Pos", "Ref", "Alt", "Qual", sample_ids)

  # the details is needed to store the genomic coordinates of the variants
  details <- genotype_data[, c("Chrom", "Pos", "Ref", "Alt", "Qual")]
  snps <- as.matrix(subset(genotype_data, select = -(1L:5L)))
  snps[snps == "0/0"]                 <- "0"
  snps[snps == "1/1"]                 <- "1"
  snps[snps == "0/1" | snps == "1/0"] <- "2"
  snps[snps == "./." | snps == ".|."] <- NA
  snps <- apply(snps, 2L, as.integer)
  meta <- add_metadata(sample_ids, meta_file)
  meta[["percentage_missing_sites"]] <- colSums(is.na(snps)) / nrow(snps)
  details[["percentage_missing_samples"]] <- rowSums(is.na(snps)) / ncol(snps)

  # adding the annotation data to the details table to associate each SNPs to
  # its gene of origin together with that gene's function.
  details[["gene"]] <- get_gene_annotation(
    genomic_coordinates = details[, c("Chrom", "Pos")],
    go = go,
    bed = bed,
    num_cores = num_threads
  )

  # we created the SNPdata class to handle easily the combined set of all the
  # data needed for downstream analyses.
  snp_table <- list(
    meta = meta,
    details = details,
    GT = snps,
    vcf = vcf_file
  )
  class(snp_table) <- "SNPdata"
  snp_table
}
