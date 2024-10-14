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
