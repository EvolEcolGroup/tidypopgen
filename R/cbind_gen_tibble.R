#' Combine a gen_tibble to a data.frame or tibble by column
#'
#' A [cbind()] method to merge `gen_tibble` objects with data.frames and normal
#' tibbles. Whilst this works, it is not ideal as it does not check the order of
#' the tables, and we suggest that you use [`dplyr::left_join()`] instead. Note
#' that `cbind` will not combine two `gen_tibbles` (i.e. it will NOT combine
#' markers for the same individuals)
#'
#' @param ... a gen_tibble and a data.frame or tibble
#' @param deparse.level an integer controlling the construction of column names.
#'   See \code{\link[base]{cbind}} for details.
#' @returns a `gen_tibble`
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # Create a dataframe to combine with the gen_tibble
#' df <- data.frame(region = c("A", "A", "B", "B", "A", "B", "B"))
#'
#' # Combine the gen_tibble with the dataframe
#' example_gt <- cbind(example_gt, df)
cbind.gen_tbl <- function(..., deparse.level = 1) { # nolint
  # send it to the data.frame method
  out <- cbind.data.frame(..., deparse.level = deparse.level)
  # sort out the lost class information
  class(out) <- c("gen_tbl", "tbl_df", "tbl", "data.frame")
  return(out)
}
