#' Principal Component Analysis for `gen_tibble` objects
#'
#' There are a number of PCA methods available for `gen_tibble` objects. They
#' are mostly designed to work on very large datasets, so they only compute a
#' limited number of components. For smaller datasets, `gt_partialSVD` allows
#' the use of partial (truncated) SVD to fit the PCA; this method is suitable
#' when the number of individuals is much smaller than the number of loci. For
#' larger dataset, `gt_randomSVD` is more appropriate. Finally, there is a
#' method specifically designed for dealing with LD in large datasets,
#' `gt_autoSVD`. Whilst this is arguably the best option, it is somewhat data
#' hungry, and so only suitable for very large datasets (hundreds of individuals
#' with several hundred thousands markers, or larger).
#'
#' NOTE: using gt_pca_autoSVD with a small dataset will likely cause an error,
#' see man page for details.
#'
#' NOTE: monomorphic markers must be removed before PCA is computed. The error
#' message 'Error: some variables have zero scaling; remove them before
#' attempting to scale.' indicates that monomorphic markers are present.
#'
#' @name gt_pca
NULL


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

square_frobenious <- function(
  X,
  ind.row = bigstatsr::rows_along(X), # nolint
  ind.col = bigstatsr::cols_along(X), # nolint
  center = rep(0, length(ind.col)),
  scale = rep(1, length(ind.col))
) {
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
#'
