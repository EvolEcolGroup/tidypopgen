#' PCA for `gen_tibble` objects by randomized partial SVD
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is deprecated: use [gt_pca()] with `method = "random"`, which
#' runs the same algorithm (and needs complete or imputed genotypes).
#'
#' This function performs Principal Component Analysis on a `gen_tibble`,
#' by randomised partial SVD based on the
#' algorithm in RSpectra (by Yixuan Qiu and Jiali Mei).\cr
#' This algorithm is linear in time in all dimensions and is very
#' memory-efficient. Thus, it can be used on very large big.matrices.
#' This function is a wrapper
#' for [bigstatsr::big_randomSVD()]
#'
#' With up to a few thousand individuals, [gt_pca_partialSVD()] is usually
#' much faster: on 1,043 individuals and 643,733 loci, this function took
#' about 175 s on one core (with the default `tol`), against about 31 s for
#' [gt_pca_partialSVD()]. This algorithm is most useful when there are too
#' many individuals for [gt_pca_partialSVD()], whose memory grows with the
#' square of the number of individuals. Using several cores (`n_cores`), or a
#' larger `tol`, also reduces its time.
#'
#' @param x a `gen_tibble` object
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
#'   and compute the product of `X` (subsetted and scaled) with `x`.
#' @param fun_cprod Same as `fun.prod`, but for the *transpose* of `X`.
#' @returns a `gt_pca` object, which is a subclass of `bigSVD`; this is
#' an S3 list with elements:
#' A named list (an S3 class "big_SVD") of
#' - `d`, the eigenvalues (singular values, i.e. as variances),
#' - `u`, the scores for each sample on each component (the left singular
#'    vectors)
#' - `v`, the loadings (the right singular vectors)
#' - `center`, the centering vector,
#' - `scale`, the scaling vector,
#' - `method`, a string defining the method (in this case 'randomSVD'),
#' - `call`, the call that generated the object.
#'
#' @note Rather than accessing these elements directly, it is better to use
#' `tidy` and `augment`. See [`gt_pca_tidiers`].
#'
#' @note Monomorphic markers must be removed before PCA is computed. The error
#' message 'Error: some variables have zero scaling; remove them before
#' attempting to scale.' indicates that monomorphic markers are present.
#' @export
#' @seealso [bigstatsr::big_randomSVD()] which this function wraps.
#' @examplesIf all(rlang::is_installed(c("RhpcBLASctl", "data.table")))
#' \dontshow{
#' data.table::setDTthreads(2)
#' RhpcBLASctl::blas_set_num_threads(2)
#' RhpcBLASctl::omp_set_num_threads(2)
#' }
#' vcf_path <-
#'   system.file("extdata", "anolis",
#'     "punctatus_t70_s10_n46_filtered.recode.vcf.gz",
#'     package = "tidypopgen"
#'   )
#' anole_gt <-
#'   gen_tibble(vcf_path, quiet = TRUE, backingfile = tempfile("anolis_"))
#'
#' # Remove monomorphic loci and impute
#' anole_gt <- anole_gt %>% select_loci_if(loci_maf(genotypes) > 0)
#' anole_gt <- gt_impute_simple(anole_gt, method = "mode")
#'
#' # Create PCA object, including total variance
#' gt_pca_randomSVD(anole_gt, k = 10, total_var = TRUE)
#' @family gt_pca_functions
# nolint start
gt_pca_randomSVD <- function(
  # nolint end
  x,
  k = 10,
  fun_scaling = bigsnpr::snp_scaleBinom(),
  tol = 1e-4,
  verbose = FALSE,
  n_cores = 1,
  fun_prod = bigstatsr::big_prodVec,
  fun_cprod = bigstatsr::big_cprodVec,
  total_var = TRUE
) {
  lifecycle::deprecate_soft(
    "1.1.0", "gt_pca_randomSVD()", "gt_pca(method = 'random')"
  )
  this_svd <- pca_random_svd(x,
    k = k, fun_scaling = fun_scaling, tol = tol, verbose = verbose,
    n_cores = n_cores, fun_prod = fun_prod, fun_cprod = fun_cprod,
    total_var = total_var
  )
  this_svd$call <- match.call()
  this_svd
}
