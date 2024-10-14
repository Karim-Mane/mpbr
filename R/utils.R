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
#'
remove_snps_from_vcf <- function(vcf, loci_to_be_retained, path, index = 1L) {
  checkmate::assert_file_exists(vcf)
  checkmate::assert_character(loci_to_be_retained, any.missing = FALSE,
                              null.ok = FALSE)
  checkmate::assert_directory_exists(path)
  checkmate::assert_numeric(index, lower = 0L, len = 1L, null.ok = FALSE,
                            any.missing = FALSE)
  target_loci <- file.path(path, loci_to_be_retained)
  header <- file.path(path, "Header.txt")
  body <- file.path(path, "Body.txt")
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

  as.character(filtered_vcf)
}
