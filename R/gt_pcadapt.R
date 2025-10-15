#' pcadapt analysis on a `gen_tibble` object
#'
#' pcadapt is an algorithm that detects genetic markers under selection. It is
#' based on the principal component analysis (PCA) of the genotypes of the
#' individuals. The method is described in Luu et al. (2017). See the R package
#' `pcadapt`, which provides extensive documentation and examples.
#'
#' Internally, this function uses the `snp_pcadapt` function from the `bigsnpr`
#' package.
#' @references Luu, K., Bazin, E., Blum, M. G. B., & François, O. (2017).
#'   pcadapt: an R package for genome scans for selection based on principal
#'   component analysis. Molecular Ecology Resources, 17(1), 67–77.
#' @param x A `gen_tibble` object.
#' @param pca a [`gt_pca`] object, as returned by `gt_pca_partialSVD()` or
#'   `gt_pca_randomSVD()`.
#' @param k Number of principal components to use in the analysis.
#' @param n_cores Number of cores to use.
#' @returns An object of subclass `gt_pcadapt`, a subclass of `mhtest`.
#' @export
#' @seealso [bigsnpr::snp_pcadapt()] which this function wraps.
#' @examplesIf all(rlang::is_installed(c("RhpcBLASctl", "data.table")))
#' \dontshow{
#' data.table::setDTthreads(2)
#' RhpcBLASctl::blas_set_num_threads(2)
#' RhpcBLASctl::omp_set_num_threads(2)
#' }
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
#' # Create PCA object
#' pca <- gt_pca_partialSVD(lobsters)
#'
#' # Create a gt_pcadapt object
#' gt_pcadapt(lobsters, pca, k = 2)
gt_pcadapt <- function(x, pca, k, n_cores = 1) {
  stopifnot_gen_tibble(x)
  if (n_cores > 1) {
    # Remove checking for two levels of parallelism
    .old_opt <- getOption("bigstatsr.check.parallel.blas", TRUE)
    options(bigstatsr.check.parallel.blas = FALSE)
    on.exit(options(bigstatsr.check.parallel.blas = .old_opt), add = TRUE)
  }
  if (!inherits(pca, "gt_pca")) {
    stop("pca must be a gt_pca object")
  }
  # check that k is a scalar
  if (!is.numeric(k) || length(k) != 1) {
    stop("k must be a scalar")
  }
  # check that k is not larger than the number of components in pca
  if (k > ncol(pca$v)) {
    stop("K must be less than or equal to the number of components in pca")
  }
  # set imputation if needed
  if (gt_has_imputed(x) && gt_uses_imputed(x) == FALSE) {
    gt_set_imputed(x, set = TRUE)
    on.exit(gt_set_imputed(x, set = FALSE), add = TRUE)
  }

  # Run the analysis
  res <- bigsnpr::snp_pcadapt(
    G = .gt_get_fbm(x),
    U.row = pca$u[, 1:k, drop = FALSE],
    ind.row = .gt_fbm_rows(x),
    ind.col = .gt_fbm_cols(x),
    ncores = n_cores
  )

  # add the loci table to the object
  attr(res, "loci") <- show_loci(x)
  class(res) <- c("gt_pcadapt", class(res))
  # Return the result
  return(res)
}
