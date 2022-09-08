#' Helper function
#'
#' @description On-the-fly NULLing of "none" or "" entries.
n2n <- function(x) {
  if (x == "none" || x == "") {
    NULL
  } else {
    x
  }
}
