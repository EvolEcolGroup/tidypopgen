#' PCA controlling for LD for `gen_tibble` objects
#'
#' This function performs Principal Component Analysis on a `gen_tibble`,
#' using a fast truncated SVD with initial pruning and then iterative removal
#' of long-range LD regions. This function is a wrapper for [bigsnpr::snp_autoSVD()]
#'
#' @param x a `gen_tbl` object
#' @param k Number of singular vectors/values to compute. Default is `10`.
#'   **This algorithm should be used to compute a few singular vectors/values.**
#' @param fun_scaling Usually this  can be left unset, as it defaults to
#' [bigsnpr::snp_scaleBinom()], which is the appropriate function for biallelic SNPs.
#' Alternatively it is possible to use  custom function
#' (see [bigsnpr::snp_autoSVD()] for details.
#' @param roll_size Radius of rolling windows to smooth log-p-values.
#'   Default is `50`.
#' @param int_min_size Minimum number of consecutive outlier SNPs
#'   in order to be reported as long-range LD region. Default is `20`.
#' @param verbose Output some information on the iterations? Default is `TRUE`.
#' @param thr_r2 Threshold over the squared correlation between two SNPs.
#'   Default is `0.2`. Use `NA` if you want to skip the clumping step.
#' size
#' @param use_positions a boolean on whether the position is used to define `size`,
#' or whether the size should be in number of SNPs. Default is TRUE
#' @param size For one SNP, window size around this SNP to compute correlations.
#'  Default is 100 / thr_r2 for clumping (0.2 -> 500; 0.1 -> 1000; 0.5 -> 200).
#'  If not providing infos.pos (NULL, the default), this is a window in number
#'  of SNPs, otherwise it is a window in kb (genetic distance). I recommend
#'  that you provide the positions if available.
#' @param alpha_tukey Default is `0.1`. The type-I error rate in outlier
#'   detection (that is further corrected for multiple testing).
#' @param min_mac Minimum minor allele count (MAC) for variants to be included.
#'   Default is `10`.
#' @param max_iter Maximum number of iterations of outlier detection.
#'   Default is `5`.
#' @param ncores Number of cores used. Default doesn't use parallelism.
#' You may use [bigstatsr::nb_cores()].
#' @returns a `gt_pca` object, which is a subclass of `bigSVD`; this is
#' an S3 list with elements:
#' A named list (an S3 class "big_SVD") of
#' - `d`, the eigenvalues (singular values, i.e. as variances),
#' - `u`, the scores for each sample on each component (the left singular vectors)
#' - `v`, the loadings (the right singular vectors)
#' - `center`, the centering vector,
#' - `scale`, the scaling vector,
#' - `method`, a string defining the method (in this case 'autoSVD'),
#' - `call`, the call that generated the object.
#'
#' Note: rather than accessing these elements directly, it is better to use
#' `tidy` and `augment`. See [`gt_pca_tidiers`].
#'
#' @export

## Look at to manipulate ellipses when passing arguments
#https://stackoverflow.com/questions/60338114/updating-values-of-three-dot-ellipsis-in-r

gt_pca_autoSVD <- function(x, k = 10,
                           fun_scaling = bigsnpr::snp_scaleBinom(),
                           thr_r2 = 0.2,
                           use_positions = TRUE,
                           size = 100 / thr_r2,
                           roll_size = 50,
                           int_min_size = 20,
                           alpha_tukey = 0.05,
                           min_mac = 10,
                           max_iter = 5,
                           ncores = 1,
                           verbose = TRUE) {

  X <- attr(x$genotypes,"bigsnp") # convenient pointer
  x_ind_col <- show_loci(x)$big_index
  x_ind_row <- vctrs::vec_data(x$genotypes)
  if (is.character(show_loci(x)$chromosome)){
    infos_chr <- as.numeric(factor(show_loci(x)$chromosome))
  } else {
    infos_chr <- show_loci(x)$chromosome
  }
  this_svd  <- bigsnpr::snp_autoSVD(X$genotypes,
                      infos.chr = infos_chr,
                      infos.pos = if (use_positions) {show_loci(x)$position} else {NULL},
                      ind.row = vctrs::vec_data(x$genotypes),
                      ind.col = show_loci(x)$big_index,
                      fun.scaling = fun_scaling,
                      thr.r2 = thr_r2,
                      size = size,
                      k = k,
                      roll.size = roll_size,
                      int.min.size = int_min_size,
                      alpha.tukey = alpha_tukey,
                      min.mac = min_mac,
                      max.iter = max_iter,
                      ncores = ncores,
                      verbose = verbose)
  # add names to the scores (to match them to data later)
  rownames(this_svd$u)<-x$id
  this_svd$method <- "autoSVD"
  this_svd$call <- match.call()
  class(this_svd) <- c("gt_pca", class(this_svd))
  this_svd
}
