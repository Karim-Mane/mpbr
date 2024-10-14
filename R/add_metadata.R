#' Build the sample metadata table
#'
#' @param sample_ids A vector of sample IDs. This should be in the same order
#'    as they appear in the input VCF file.
#' @param metadata The path to the sample metadata file
#'
#' @return An object of class `data.frame` that contains the sample metadata
#' @keywords internal
#'
#' @details
#' In addition to generating the sample metadata file, the function will also
#' create a file, named as `samples_in_metadata_not_in_vcf.txt`, if there are
#' samples in the metadata file that were not found in the VCF file.
#'
add_metadata <- function(sample_ids, metadata) {
  checkmate::assert_vector(sample_ids, min.len = 1L, null.ok = FALSE)
  checkmate::assert_file_exists(metadata)
  if (grepl(".RDS", basename(metadata))) {
    meta <- readRDS(metadata)
  } else {
    meta <- data.table::fread(metadata, key = "sample", nThread = 4L)
  }
  samples <- meta[["sample"]]

  # sample from the VCF file must match with those in the metadata file
  are_in_meta_file <- sample_ids %in% samples
  if (!all(are_in_meta_file)) {
    warning("Incomplete metadata - the following samples in the VCF",
            " file are not found in metadata file:\n",
            toString(sample_ids[!are_in_meta_file]),
            call. = FALSE)
  }

  # samples from the metadata file should also match with the ones from
  # the input VCF file
  are_in_vcf_file <- samples %in% sample_ids
  if (!all(are_in_vcf_file)) {
    warning(
      "The following samples are removed from metadata file as ",
      "they are not found in the VCF file: ",
      toString(samples[!are_in_vcf_file]),
      call. = FALSE
    )
    meta <- meta[-(!are_in_vcf_file), ]
  }

  # joining the samples IDs from the VCF file with the sample metadata
  meta <- data.frame(sample = samples) %>%
    dplyr::left_join(meta, by = "sample")

  meta
}
