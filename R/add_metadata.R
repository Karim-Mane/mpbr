#' Build the sample metadata table
#'
#' @param sample_ids a `vector` of sample IDs. This should be in the same order
#'    as they appear in the input VCF file.
#' @param metadata the path to the sample metadata file
#' @param output_dir the path to the output directory
#'
#' @return an object of class `data.frame` that contains the sample metadata
#' @keywords internal
#' @noRd
#'
#' @details
#' In addition to generating the sample metadata file, the function will also
#' create a file, named as `samples_in_metadata_not_in_vcf.txt`, if there are
#' samples in the metadata file that were not found in the VCF file.
#'
add_metadata <- function(sample_ids, metadata, output_dir) {
  checkmate::assert_data_frame(sample_ids, min.rows = 1L, min.cols = 1L,
                               null.ok = FALSE)
  checkmate::assert_file_exists(metadata)
  checkmate::assert_directory_exists(output_dir)
  if (grepl(".RDS", basename(metadata))) {
    meta <- readRDS(metadata)
  } else {
    data.table::fread(metadata, key = "sample", nThread = 4L)
  }
  
  samples <- meta[["sample"]]
  
  # sample from the VCF file must match with those in the metadata file
  are_in_meta_file <- sample_ids[["sample"]] %in% samples
  if (!all(are_in_meta_file)) {
    warning("Incomplete metadata - the following samples in the VCF",
            "file are not found in metadata file:\n",
            glue::glue_collapse(sample_ids[["sample"]][!are_in_meta_file],
                                sep = ", "),
            call. = FALSE)
    data.table::fwrite(sample_ids[["sample"]][!are_in_meta_file],
                       file.path(output_dir,
                                 "samples_in_metadata_not_in_vcf.txt"),
                       row.names = FALSE, col.names = FALSE, quote = FALSE)
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