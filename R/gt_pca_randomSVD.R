#' PCA for `gen_tibble` objects by randomized partial SVD
#'
#' This function performs Principal Component Analysis on a `gen_tibble`,
#' by randomised partial SVD based on the
#' algorithm in RSpectra (by Yixuan Qiu and Jiali Mei).\cr
#' This algorithm is linear in time in all dimensions and is very
#' memory-efficient. Thus, it can be used on very large big.matrices.
#' This function is a wrapper
#' for [bigstatsr::big_randomSVD()]
#'
#' @param x a `gen_tbl` object
#' @param k Number of singular vectors/values to compute. Default is `10`.
#'   **This algorithm should be used to compute a few singular vectors/values.**
#' @param fun_scaling Usually this  can be left unset, as it defaults to
#' [bigsnpr::snp_scaleBinom()], which is the appropriate function for biallelic SNPs.
#' Alternatively it is possible to use  custom function
#' (see [bigsnpr::snp_autoSVD()] for details.
#' @param n_cores Number of cores used.
#' @param tol Precision parameter of [svds][RSpectra::svds]. Default is `1e-4`.
#' @param verbose Should some progress be printed? Default is `FALSE`.
#' @param fun_prod Function that takes 6 arguments (in this order):
#'  - a matrix-like object `X`,
#'  - a vector `x`,
#'  - a vector of row indices `ind.row` of `X`,
#'  - a vector of column indices `ind.col` of `X`,
#'  - a vector of column centers (corresponding to `ind.col`),
#'  - a vector of column scales (corresponding to `ind.col`),
#'  and compute the product of `X` (subsetted and scaled) with `x`.
#' @param fun_cprod Same as `fun.prod`, but for the *transpose* of `X`.

#' @returns a `gt_pca` object, which is a subclass of `bigSVD`; this is
#' an S3 list with elements:
#' A named list (an S3 class "big_SVD") of
#' - `d`, the eigenvalues (singular values, i.e. as variances),
#' - `u`, the scores for each sample on each component (the left singular vectors)
#' - `v`, the loadings (the right singular vectors)
#' - `center`, the centering vector,
#' - `scale`, the scaling vector,
#' - `method`, a string defining the method (in this case 'randomSVD'),
#' - `call`, the call that generated the object.
#'
#' Note: rather than accessing these elements directly, it is better to use
#' `tidy` and `augment`. See [`gt_pca_tidiers`].
#' @export


gt_pca_randomSVD <- function(x, k = 10,
                              fun_scaling = bigsnpr::snp_scaleBinom(),
                              tol = 1e-4,
                              verbose = FALSE,
                              n_cores = 1,
                              fun_prod = bigstatsr::big_prodVec,
                              fun_cprod = bigstatsr::big_cprodVec
                              ) {
  X <- attr(x$genotypes,"bigsnp") # convenient pointer
  x_ind_col <- show_loci(x)$big_index
  x_ind_row <- vctrs::vec_data(x$genotypes)
  if (is.character(show_loci(x)$chromosome)){
    infos_chr <- as.numeric(factor(show_loci(x)$chromosome))
  } else {
    infos_chr <- show_loci(x)$chromosome
  }
  this_svd  <- bigstatsr::big_randomSVD(X$genotypes,
                                  k=k,
                                    ind.row = vctrs::vec_data(x$genotypes),
                                    ind.col = show_loci(x)$big_index,
                                    fun.scaling = fun_scaling,
                                    tol = tol,
                                  verbose = verbose,
                                  n_cores = n_cores,
                                  fun.prod = fun_prod,
                                  fun.cprod = fun_cprod) # TODO check that this is correct and expose it, maybe creat convenience function to get the values
  # add names to the scores (to match them to data later)
  rownames(this_svd$u)<-x$id
  rownames(this_svd$v) <- show_loci_names(x)
  this_svd$method <- "randomSVD"
  this_svd$call <- match.call()
  class(this_svd) <- c("gt_pca", class(this_svd))
  this_svd
}
