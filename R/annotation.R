#' Associate each SNPs to a gene id, name and description based on the gene
#' they fall under
#'
#' @param target_gtf A data frame with the gene annotation
#' @param genomic_coordinates A data frame with the SNPs genomic coordinates
#'
#' @return A vector with gene annotation. This should be of the same length as
#'    the number of SNPs.
#' @keywords internal
#'
gene_annotation <- function(target_gtf, genomic_coordinates) {
  checkmate::assert_data_frame(target_gtf, min.rows = 1L, min.cols = 1L,
                               null.ok = FALSE)

  names(genomic_coordinates) <- c("chrom", "start")
  genomic_coordinates[["end"]] <- genomic_coordinates[["start"]]
  subject <- IRanges::IRanges(target_gtf[["start"]], target_gtf[["end"]])
  query <- IRanges::IRanges(
    genomic_coordinates[["start"]],
    genomic_coordinates[["end"]]
  )
  my_overlaps <- data.table::data.table(
    as.matrix(GenomicRanges::findOverlaps(
      query,
      subject,
      type = "within"))
  )
  my_overlaps[["gene_id"]] <- target_gtf[["gene_id"]][my_overlaps[["subjectHits"]]] #nolint: line_length_linter
  my_overlaps[["gene_name"]] <- target_gtf[["gene_name"]][my_overlaps[["subjectHits"]]] #nolint: line_length_linter
  my_overlaps[["gene_desc"]] <- target_gtf[["gene_desc"]][my_overlaps[["subjectHits"]]] #nolint: line_length_linter
  my_overlaps <- my_overlaps[, lapply(.SD, paste, collapse = ":"), by = queryHits]
  my_overlaps[["queryHits"]] <- NULL
  my_overlaps[["subjectHits"]] <- NULL
  
  genomic_coordinates[["end"]] <- NULL
  genomic_coordinates <- cbind(genomic_coordinates, my_overlaps)

  return(genomic_coordinates)
}

#' Clean gene names
#'
#' @param y The gene name
#'
#' @return A string with the cleaned gene name
#' @keywords internal
#'
get_clean_name <- function(y) {
  checkmate::assert_character(y, any.missing = FALSE, null.ok = FALSE)
  unlist(strsplit(
    unlist(strsplit(y, "ID=", fixed = TRUE))[[2L]],
    ";", fixed = TRUE))[[1L]]
}

get_gene_id_name_desc <- function(bed) {
  # extract gene IDs
  gene_id <- as.character(
    lapply(bed[["V10"]], function(x) {
      unlist(strsplit(x, ";", fixed = TRUE))[[1]]
    })
  )
  gene_id <- as.character(
    lapply(gene_id, function(x) {unlist(strsplit(x, "=", fixed = TRUE))[[2]]})
  )
  
  # extract gene names
  gene_name <- as.character(
    lapply(bed[["V10"]], function(x) {
      unlist(strsplit(x, ";", fixed = TRUE))[[2]]
    })
  )
  gene_name <- as.character(
    lapply(gene_name, function(x) {unlist(strsplit(x, "=", fixed = TRUE))[[2]]})
  )
  
  # extract gene descriptions
  gene_desc <- as.character(
    lapply(bed[["V10"]], function(x) {
      unlist(strsplit(x, ";", fixed = TRUE))[[3]]
    })
  )
  gene_desc <- as.character(
    lapply(gene_desc, function(x) {unlist(strsplit(x, "=", fixed = TRUE))[[2]]})
  )
  
  return(data.frame(
    gene_id = gene_id,
    gene_name = gene_name,
    gene_desc = gene_desc,
    stringsAsFactors = FALSE
  ))
}

#' Add gene ontology and names annotation details to every SNPs in the table
#' that contains their genomic coordinates
#'
#' @param genomic_coordinates The table with SNPs genomic coordinates
#' @param go A data frame with the gene ontology annotation details
#' @param bed A data frame with the gene name annotation
#' @param nthreads The number of cores to be used
#'
#' @return An object of type `data.frame` with the SNPs genomic coordinates and
#'    their corresponding annotation details
#' @keywords internal
#'
get_gene_annotation <- function(genomic_coordinates, go, bed, nthreads = 4L) {
  checkmate::assert_data_frame(genomic_coordinates, min.rows = 1L,
                               min.cols = 1L, null.ok = FALSE)
  checkmate::assert_data_frame(go, min.rows = 1L,
                               min.cols = 1L, null.ok = FALSE)
  checkmate::assert_data_frame(bed, min.rows = 1L,
                               min.cols = 1L, null.ok = FALSE)

  # define an IRanges object with the gene genomic coordinates
  # genes_genomic_coordinates <- GenomicRanges::GRanges(
  #   seqnames = bed[["V1"]],
  #   ranges = IRanges::IRanges(start = bed[["V2"]], end = bed[["V3"]])
  # )
  # genes_genomic_coordinates <- sort(genes_genomic_coordinates)
  genes_genomic_coordinates <- subset(bed, select = 1:3)
  names(genes_genomic_coordinates) <- c("chrom", "start", "end")

  # extract the gene IDs, names, and description
  gene_annot <- get_gene_id_name_desc(bed)
  gene_annot <- cbind(genes_genomic_coordinates, gene_annot)

  # select the gene ID and description columns from the gene ontology
  go <- subset(go, select = c(2, 10))
  names(go) <- c("gene_id", "gene_name")
  
  chrom <- gene_id <- gene_name <- start <- end <- NULL # nolint: object_name_linter
  # joint gene annotation and ontology
  gene_annot_onto <- gene_annot |>
    dplyr::left_join(go, by = "gene_id", relationship = "many-to-many")
  gene_annot_onto <- dplyr::distinct(gene_annot_onto, .keep_all = TRUE)
  gene_annot_onto[["gene_desc"]] <- gene_annot_onto[["gene_name.y"]]
  gene_annot_onto[["gene_name.y"]] <- NULL
  names(gene_annot_onto)[[5]] <- "gene_name"
  
  # add annotation to each SNP
  genomic_coordinates[["Pos"]] <- as.numeric(genomic_coordinates[["Pos"]])
  resultat <- gene_annotation(
    target_gtf = gene_annot_onto,
    genomic_coordinates = genomic_coordinates
  )
  resultat <- data.frame(
    apply(resultat, 2, function(x) {gsub("NA:", "", x, fixed = TRUE)})
  )
  return(resultat)
}


