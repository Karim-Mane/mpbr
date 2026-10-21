#' Print information about the <SNPdata> object
#'
#' @param x An object of class <SNPdata>
#' @param ... further arguments passed to or from other methods.
#'
#' @returns prints out the structure of the <SNPdata> object on the console
#'
#' @export
#'
print.SNPdata <- function(x, ...) {
  checkmate::assert_class(x, "SNPdata", null.ok = FALSE)
  cli::cli_alert_info("sample metadata")
  print(utils::head(x[["meta"]]))
  cli::cli_alert_info("\nSNPs genomic coordinates and annotations")
  print(utils::head(x[["details"]]))
  cli::cli_inform(c(
    "*" = "\nNumber of samples: {.val {dim(x[[\"GT\"]])[[2]]}}",
    "*" = "Number of SNPs: {.val {dim(x[[\"GT\"]])[[1]]}}",
    "*" = "Input VCF file: {.file {x[[\"vcf\"]]}}"
  ))
}