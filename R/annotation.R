#' Associate each SNPs to a gene name on which it belongs
#'
#' @param target_gtf a `data.frame` with the gene annotation
#' @param genomic_coordinates a `data.frame` with the SNPs genomic coordinates
#'
#' @return a `vector` with gene annotation. This should be of the same length as
#'    the number of SNPs.
#' @keywords internal
#'
gene_annotation <- function(target_gtf, genomic_coordinates) {
  checkmate::assert_data_frame(target_gtf, min.rows = 1L, min.cols = 1L,
                               null.ok = FALSE)
  
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
  gene <- queryHits <- NULL # nolint: object_name_linter
  the_genes <- my_overlaps %>%
    dplyr::group_by(queryHits) %>%
    dplyr::distinct(queryHits, gene) %>%
    dplyr::filter(!is.na(gene)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(queryHits) %>%
    dplyr::summarize(gene = glue::glue_collapse(gene, sep = ":"))
  genomic_coordinates[["gene"]]    <- NA
  genomic_coordinates[["gene"]][the_genes[["queryHits"]]] <- the_genes[["gene"]]
  
  genomic_coordinates[["gene"]]
}

#' clean gene names
#'
#' @param y the gene name
#'
#' @return a string with the cleaned gene name
#' @keywords internal
#'
get_clean_name <- function(y) {
  checkmate::assert_character(y, any.missing = FALSE, null.ok = FALSE)
  unlist(strsplit(
    unlist(strsplit(y, "ID=", fixed = TRUE))[[2L]],
    ";", fixed = TRUE))[[1L]]
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
  go           <- subset(go, select = c(2, 10))
  names(go)    <- c("gene_id", "gene_name")
  
  chrom <- gene_id <- gene_name <- start <- end <- NULL # nolint: object_name_linter
  test         <- genes %>% dplyr::left_join(go, by = "gene_id",
                                             relationship = "many-to-many")
  test         <- dplyr::distinct(test, chrom, start, end, gene_id, gene_name)
  resultat     <- gene_annotation(test, genomic_coordinates)
  resultat     <- gsub("NA:", "", resultat, fixed = TRUE)
  resultat
}