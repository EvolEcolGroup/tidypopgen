#' Combine a gen_tibble to a data.frame or tibble by column
#'
#' A [cbind()] method for `gen_tibble` objects. Whilst this works, it is not
#' ideal, and we suggest that you use
#' [`dplyr::bind_cols()`] instead.
#'
#' @param ... a gen_tibble and a data.frame or tibble
#' @param deparse.level an integer controlling the construction of column names.
#' See \code{\link[base]{cbind}} for details.
#' @returns a `gen_tibble`
#' @export

cbind.gen_tbl <- function(..., deparse.level = 1) { #nolint
  # send it to the next (data.frame) method
  out <- NextMethod(generic = "cbind", ..., deparse.level = deparse.level)
  # sort out the lost class information
  class(out) <- c("gen_tbl", "tbl_df", "tbl", "data.frame")
  return(out)
}
