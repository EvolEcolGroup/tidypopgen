#' Compute the square of the Frobenius norm of a matrix
#'
#' This function computes the square of the Frobenius norm of a matrix, which is
#' the sum of the squares of the matrix elements, which provides a measure of
#' of the total variance of the matrix.
#' The code used here was outlined in an issue of bigstatsr by @privefl:
#' https://github.com/privefl/bigstatsr/issues/83
#'
#' @param X An FBM matrix
#' @param ind.row A vector of row indices
#' @param ind.col A vector of column indices
#' @param center A vector of centering values (i.e the means of the genotype
#'   counts)
#' @param scale A vector of scaling values (i.e the standard deviations of the
#'   genotype counts)
#' @return The square of the Frobenius norm of the matrix
#' @keywords internal
#' @noRd
square_frobenius <- function(
    X, # nolint start
    ind.row = bigstatsr::rows_along(X),
    ind.col = bigstatsr::cols_along(X), # nolint end
    center = rep(0, length(ind.col)),
    scale = rep(1, length(ind.col))) {
  if (length(center) != length(ind.col) || length(scale) != length(ind.col)) {
    stop(paste(
      "center and scale must be the same length as the number of",
      "columns in the matrix"
    ))
  }
  stats <- bigstatsr::big_colstats(X, ind.row, ind.col)
  n <- length(ind.row)
  sum(((n - 1) * stats$var + n * (stats$sum / n - center)^2) / scale^2)
}
