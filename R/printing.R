#' Print the `SNPdata` object
#'
#' @param x An object of class `SNPdata`
#' @param ... further arguments passed to or from other methods.
#'
#' @returns prints the structure of the `SNPdata` object
#'
#' @export
#'
print.SNPdata <- function(x, ...) {
  checkmate::assert_class(x, "SNPdata", null.ok = FALSE)
  message("metadata...")
  print(utils::head(x[["meta"]]))
  message("\ndetails on genomic coordinates...")
  print(utils::head(x[["details"]]))
  message(sprintf("\nData contains: %d samples for %d snp loci",
                  dim(x[["GT"]])[[2L]],
                  dim(x[["GT"]])[[1L]]))
  message(sprintf("\nData is generated from: %s", x[["vcf"]]))
  
  # make use of the following:
  # cli::cat_line("This is ", "a ", "line of text.", col = "red")
  # cli::cat_bullet(letters[1:5])
  # cli::cat_bullet(letters[1:5], bullet = "tick", bullet_col = "green")
  # cli::cat_rule()
}