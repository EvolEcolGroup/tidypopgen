#' Combine a gen_tibble to a data.frame or tibble by column
#'
#' A [cbind()] method to merge `gen_tibble` objects with data.frames and normal
#' tibbles. Whilst this works, it is not ideal, and we suggest that you use
#' [`dplyr::bind_cols()`] instead. Note that `cbind` will not combine two
#' `gen_tibbles` (i.e. it will NOT combine markers for the same individuals)
#'
#' @param ... a gen_tibble and a data.frame or tibble
#' @param deparse.level an integer controlling the construction of column names.
#'   See \code{\link[base]{cbind}} for details.
#' @returns a `gen_tibble`
#' @export

cbind.gen_tbl <- function(..., deparse.level = 1) { # nolint
  # send it to the next (data.frame) method
  out <- NextMethod(generic = "cbind", ..., deparse.level = deparse.level)
  # sort out the lost class information
  class(out) <- c("gen_tbl", "tbl_df", "tbl", "data.frame")
  return(out)
}
