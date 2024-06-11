#' PCA for `gen_tibble` objects by partial SVD
#'
#' This function performs Principal Component Analysis on a `gen_tibble`,
#' by partial SVD through the eigen decomposition of the covariance. It works well
#' if the number of individuals is much smaller than the number of loci; otherwise,
#' [gt_pca_randomSVD()] is a better option. This function is a wrapper
#' for [bigstatsr::big_SVD()]
#'
#' @param x a `gen_tbl` object
#' @param k Number of singular vectors/values to compute. Default is `10`.
#'   **This algorithm should be used to compute a few singular vectors/values.**
#' @param fun_scaling Usually this  can be left unset, as it defaults to
#' [bigsnpr::snp_scaleBinom()], which is the appropriate function for biallelic SNPs.
#' Alternatively it is possible to use  custom function
#' (see [bigsnpr::snp_autoSVD()] for details.
#' @returns a `gt_pca` object, which is a subclass of `bigSVD`; this is
#' an S3 list with elements:
#' A named list (an S3 class "big_SVD") of
#' - `d`, the eigenvalues (singular values, i.e. as variances),
#' - `u`, the scores for each sample on each component (the left singular vectors)
#' - `v`, the loadings (the right singular vectors)
#' - `center`, the centering vector,
#' - `scale`, the scaling vector,
#' - `method`, a string defining the method (in this case 'partialSVD'),
#' - `call`, the call that generated the object.
#'
#' Note: rather than accessing these elements directly, it is better to use
#' `tidy` and `augment`. See [`gt_pca_tidiers`].
#' @export


gt_pca_partialSVD <- function(x, k = 10, fun_scaling = bigsnpr::snp_scaleBinom()
                              ) {
  if (gt_has_imputed(x) && gt_uses_imputed(x)==FALSE){
    gt_set_imputed(x, set = TRUE)
    on.exit(gt_set_imputed(x, set = FALSE))
  }
  X <- attr(x$genotypes,"bigsnp") # convenient pointer
  x_ind_col <- show_loci(x)$big_index
  x_ind_row <- vctrs::vec_data(x$genotypes)
  if (is.character(show_loci(x)$chromosome)){
    infos_chr <- as.numeric(factor(show_loci(x)$chromosome))
  } else {
    infos_chr <- show_loci(x)$chromosome
  }
  this_svd  <- bigstatsr::big_SVD(X$genotypes,
                                  k=k,
                                    ind.row = vctrs::vec_data(x$genotypes),
                                    ind.col = show_loci(x)$big_index,
                                    fun.scaling = fun_scaling,
                                    block.size= bigstatsr::block_size(nrow(X$genotypes))) # TODO check that this is correct and expose it, maybe creat convenience function to get the values
  # add names to the scores (to match them to data later)
  rownames(this_svd$u)<-x$id
  rownames(this_svd$v) <- loci_names(x)
  this_svd$method <- "partialSVD"
  this_svd$call <- match.call()
  this_svd$loci <- show_loci(x)
  class(this_svd) <- c("gt_pca", class(this_svd))
  this_svd
}
