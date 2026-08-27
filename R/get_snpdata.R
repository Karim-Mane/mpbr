#' Create A `SNPdata` object
#'
#' The function generates the input object type require to use functions in
#' the {mpbr} package.
#'
#' @param vcf_file A character with the path to the input VCF file (required)
#' @param meta_file A character with the path to the sample metadata file
#'    (required)
#' @param output_dir A character with the path to the folder where the output
#'    and temporary files will be stored (optional)
#' @param gaf A character with the path to the file with the gene ontology
#'    annotation file (required). This can be downloaded from
#'    [PlasmoDB](https://plasmodb.org/plasmo/app/downloads/Current_Release/Pfalciparum3D7/gaf/) # nolint: line_length_linter
#' @param gff A character with the path to the gene annotation file. This file
#'    can be downloaded from
#'    [PlasmoDB](https://plasmodb.org/plasmo/app/downloads/Current_Release/Pfalciparum3D7/gff/) # nolint: line_length_linter
#' @param num_threads The number of threads to be during data processing
#'
#' @return An object of class <SNPdata> that contains the following 4 elements:
#' \enumerate{
#'   \item meta: A `data.frame` that with the sample's metadata
#'   \item details: A `data.frame` with the SNPs genomic coordinates, and their
#'      associated annotations from both the VCF and input annotation files
#'   \item GT: An integer `matrix` with the sample's genotype data where:
#'      0 represents the reference allele, 1 the alternative allele, 2 is a
#'      mixed allele, and NA a missing genotype.
#'   \item vcf: the full path to the input VCF file from which the data is
#'      generated.
#'   }
#'
#' @export
#' @examples
#' \dontrun{
#'   snpdata <- get_snpdata(
#'     vcf_file = system.file("extdata", "test_data.vcf.gz", package = "mpbr"),
#'     meta_file = system.file("extdata", "sample_metadata.RDS",
#'                             package = "mpbr"),
#'     output_dir = tempdir(),
#'     gaf = file.path("data", "PlasmoDB-68_Pfalciparum3D7_Curated_GO.gaf.gz"),
#'     gff = file.path("data", "PlasmoDB-68_Pfalciparum3D7.gff")
#'   )
#' }
get_snpdata <- function(vcf_file,
                        meta_file,
                        output_dir,
                        gaf,
                        gff,
                        num_threads = 4L) {
  checkmate::assert_file_exists(vcf_file)
  checkmate::assert_file_exists(meta_file)
  checkmate::assert_directory_exists(output_dir)
  checkmate::assert_character(gaf, any.missing = FALSE, len = 1L,
                              null.ok = TRUE)
  checkmate::assert_numeric(num_threads, lower = 1L,
                            any.missing = FALSE, null.ok = FALSE, len = 1L)
  cli::cli_progress_step("Building the genotype matrix")
  # the sample IDs will be used to create the sample metadata file.
  sample_ids <- get_sample_ids(vcf_file)

  # import the VCF file
  vcf_data <- data.table::fread(
    file = vcf_file,
    skip = "#CHROM",
    header = TRUE,
    sep = "\t",
    nThread = num_threads
  )
  # rename #CHROM to CHROM
  data.table::setnames(vcf_data, "#CHROM", "CHROM")

  # get the genotype data to be used to create the genotype matrix (GT) and the
  # details table
  genotype_data <- extract_genotype(
    vcf_data = vcf_data,
    which = "GT",
    nthreads = num_threads
  )
  names(genotype_data)[1:4] <- c("Chrom", "Pos", "Ref", "Alt")

  # convert the GT matrix into numeric
  snps <- as.matrix(subset(genotype_data, select = -(1L:4L)))
  snps[snps == "0/0" | snps == "0|0"] <- "0"
  snps[snps == "1/1" | snps == "1|1"] <- "1"
  snps[snps == "0/1" | snps == "0|1" | snps == "1|0" | snps == "1/0"] <- "2"
  snps[snps == "./." | snps == ".|."] <- NA
  snps <- apply(snps, 2L, as.integer)

  cli::cli_progress_step("Building the sample metadata table")
  meta <- add_metadata(sample_ids, meta_file)
  meta[["percentage_missing_sites"]] <- round(
    colSums(is.na(snps)) / nrow(snps),
    digits = 5
  )

  cli::cli_progress_step("Building the SNPs genomic coordinates table")
  details <- genotype_data[, c("Chrom", "Pos", "Ref", "Alt")]
  details[["percentage_missing_samples"]] <- round(
    rowSums(is.na(snps)) / ncol(snps),
    digits = 5
  )

  # import the gene ontology annotation file
  cli::cli_progress_step("Reading annotation files")
  go <- data.table::fread(
    cmd = sprintf("pigz -dc < %s", gaf),
    nThread = num_threads,
    sep = "\t",
    skip = 1
  )

  # import the GFF file
  # convert it into BED format
  # then create a gene-level only features BED
  # This is converted into BED format because the {GenomicRange}
  # package requires that file type.
  bed_file <- file.path(output_dir, "gene_annotation.bed")
  system(sprintf("gff2bed < %s > %s", gff, bed_file))
  gene_level_bed <- file.path(output_dir, "genes_only.bed")
  system(sprintf("awk '$8 ~ /_gene$/' %s > %s", bed_file, gene_level_bed))
  bed <- data.table::fread(gene_level_bed, nThread = num_threads, sep = "\t")

  # add the annotation data to the details table to associate each SNPs to
  # its gene of origin together with that gene's function
  cli::cli_progress_step("Annotating the SNPs genomic coordinates table")
  # extract SNP effect annotation and add them to the details table
  snp_effect_annot <- extract_annotation(
    vcf_data = vcf_data,
    nthreads = num_threads
  )
  names(snp_effect_annot)[1:4] <- c("Chrom", "Pos", "Ref", "Alt")
  details <- details |>
    dplyr::left_join(snp_effect_annot, by = c("Chrom", "Pos", "Ref", "Alt"))
  
  # add gene IDs, names, and description of genes under which the SNPs fall
  snps_gene_annot <- get_gene_annotation(
    genomic_coordinates = details[, c("Chrom", "Pos")],
    go = go,
    bed = bed,
    nthreads = num_threads
  )
  names(snps_gene_annot)[1:2] <- c("Chrom", "Pos")
  snps_gene_annot[["Pos"]] <- as.numeric(snps_gene_annot[["Pos"]])
  details <- details |>
    dplyr::left_join(snps_gene_annot, by = c("Chrom", "Pos"))
  
  # we created the SNPdata class to handle easily the combined set of all the
  # data needed for downstream analyses.
  cli::cli_progress_step("Building the {.cls SNPdata} object")
  snp_table <- list(
    meta = meta,
    details = details,
    GT = snps,
    vcf = vcf_file
  )
  class(snp_table) <- "SNPdata"
  return(snp_table)
}
