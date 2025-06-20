#' PCA for `gen_tibble` objects by partial SVD
#'
#' This function performs Principal Component Analysis on a `gen_tibble`, by
#' partial SVD through the eigen decomposition of the covariance. It works well
#' if the number of individuals is much smaller than the number of loci;
#' otherwise, [gt_pca_randomSVD()] is a better option. This function is a
#' wrapper for [bigstatsr::big_SVD()].
#'
#' @param x a `gen_tbl` object
#' @param k Number of singular vectors/values to compute. Default is `10`.
#'   **This algorithm should be used to compute a few singular vectors/values.**
#' @param fun_scaling Usually this  can be left unset, as it defaults to
#'   [bigsnpr::snp_scaleBinom()], which is the appropriate function for
#'   biallelic SNPs. Alternatively it is possible to use  custom function (see
#'   [bigsnpr::snp_autoSVD()] for details.
#' @param total_var a boolean indicating whether to compute the total variance
#'   of the matrix. Default is `TRUE`. Using `FALSE` will speed up computation,
#'   but the total variance will not be stored in the output (and thus it will
#'   not be possible to assign a proportion of variance explained to the
#'   components).
#' @returns a `gt_pca` object, which is a subclass of `bigSVD`; this is an S3
#'   list with elements: A named list (an S3 class "big_SVD") of
#' - `d`, the eigenvalues (singular values, i.e. as variances),
#' - `u`, the scores for each sample on each component (the left singular
#'    vectors)
#' - `v`, the loadings (the right singular vectors)
#' - `center`, the centering vector,
#' - `scale`, the scaling vector,
#' - `method`, a string defining the method (in this case 'partialSVD'),
#' - `call`, the call that generated the object.
#' - `square_frobenius`, used to compute the proportion of variance explained
#'    by the components (optional)
#'
#'   Note: rather than accessing these elements directly, it is better to use
#'   `tidy` and `augment`. See [`gt_pca_tidiers`].
#' @export
#' @examples
#' # Create a gen_tibble of lobster genotypes
#' bed_file <-
#'   system.file("extdata", "lobster", "lobster.bed", package = "tidypopgen")
#' lobsters <- gen_tibble(bed_file,
#'   backingfile = tempfile("lobsters"),
#'   quiet = TRUE
#' )
#'
#' # Remove monomorphic loci and impute
#' lobsters <- lobsters %>% select_loci_if(loci_maf(genotypes) > 0)
#' lobsters <- gt_impute_simple(lobsters, method = "mode")
#'
#' # Create PCA object, including total variance
#' gt_pca_partialSVD(lobsters,
#'   k = 10,
#'   total_var = TRUE
#' )
#' # Change number of components and exclude total variance
#' gt_pca_partialSVD(lobsters,
#'   k = 5,
#'   total_var = FALSE
#' )
# nolint start
gt_pca_partialSVD <- function(
    # nolint end
    x,
    k = 10,
    fun_scaling = bigsnpr::snp_scaleBinom(),
    total_var = TRUE) {
  if (gt_has_imputed(x) && gt_uses_imputed(x) == FALSE) {
    gt_set_imputed(x, set = TRUE)
    on.exit(gt_set_imputed(x, set = FALSE))
  }
  X <- attr(x$genotypes, "bigsnp") # convenient pointer #nolint
  x_ind_col <- show_loci(x)$big_index
  x_ind_row <- vctrs::vec_data(x$genotypes)

  this_svd <- bigstatsr::big_SVD(
    X$genotypes, # nolint
    k = k,
    ind.row = x_ind_row,
    ind.col = x_ind_col,
    fun.scaling = fun_scaling,
    block.size = bigstatsr::block_size(nrow(X$genotypes))
  ) # TODO check that this is correct and expose it, maybe create convenience
  # function to get the values
  # add names to the scores (to match them to data later)
  rownames(this_svd$u) <- x$id
  rownames(this_svd$v) <- loci_names(x)
  this_svd$method <- "partialSVD"
  this_svd$call <- match.call()
  this_svd$loci <- show_loci(x)
  class(this_svd) <- c("gt_pca", class(this_svd))
  if (total_var) {
    this_svd$square_frobenius <- square_frobenius(
      X$genotypes,
      x_ind_row,
      x_ind_col,
      center = this_svd$center,
      scale = this_svd$scale
    )
  }
  this_svd
}
