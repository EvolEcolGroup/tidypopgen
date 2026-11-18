#' pcadapt analysis on a `gen_tibble` object
#'
#' pcadapt is an algorithm that detects genetic markers under selection. It is
#' based on the principal component analysis (PCA) of the genotypes of the
#' individuals. The method is described in Luu et al. (2017). See the R package
#' `pcadapt`, which provides extensive documentation and examples.
#'
#' Internally, this function uses the `snp_pcadapt` function from the `bigsnpr`
#' package.
#'
#' The genotypes of `x` must be complete or imputed (see [gt_impute_simple()];
#' imputed genotypes are used automatically). Missing genotypes give an
#' error; as in `bigsnpr`, they are looked for in a sample of the genotypes
#' (all the individuals, at evenly spaced loci), so a few missing genotypes
#' can go unnoticed and change the results.
#' @references Luu, K., Bazin, E., Blum, M. G. B., & François, O. (2017).
#'   pcadapt: an R package for genome scans for selection based on principal
#'   component analysis. Molecular Ecology Resources, 17(1), 67–77.
#' @param x A `gen_tibble` object.
#' @param pca a [`gt_pca`] object, as returned by [gt_pca()].
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
#' # Remove monomorphic loci and impute (pcadapt needs complete genotypes)
#' lobsters <- lobsters %>% select_loci_if(loci_maf(genotypes) > 0)
#' lobsters <- gt_impute_simple(lobsters, method = "mode")
#'
#' # Create PCA object
#' pca <- gt_pca(lobsters)
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
  # check that the scores are in the same order as the individuals in x
  if (
    !is.null(rownames(pca$u)) &&
      !identical(rownames(pca$u), as.character(x$id))
  ) {
    stop("the individuals in 'x' do not correspond to those in 'pca'")
  }
  # set imputation if needed
  if (gt_has_imputed(x) && gt_uses_imputed(x) == FALSE) {
    gt_set_imputed(x, set = TRUE)
    on.exit(gt_set_imputed(x, set = FALSE), add = TRUE)
  }
  # pcadapt needs complete genotypes: bigsnpr's check samples the whole
  # backing matrix (including individuals and loci that x does not use), so
  # turn it off, and check a sample of the genotypes of x instead
  local_no_bigstatsr_checks()
  stopifnot_no_missing(x, "gt_pcadapt()", n_cores = n_cores)

  x_ind_col <- show_loci(x)$big_index
  x_ind_row <- vctrs::vec_data(x$genotypes)

  # Run the analysis
  res <- tryCatch(
    bigsnpr::snp_pcadapt(
      G = .gt_get_fbm(x),
      U.row = pca$u[, 1:k, drop = FALSE],
      ind.row = x_ind_row,
      ind.col = x_ind_col,
      ncores = n_cores
    ),
    error = function(e) {
      # loci with the same genotype in every individual have no variance,
      # and so undefined statistics, which snp_pcadapt() reports as "You
      # can't have missing values in 'U'" (#179): look for them (only after
      # this error) and say which they are
      if (grepl("missing values in 'U'", conditionMessage(e), fixed = TRUE)) {
        stopifnot_variable_loci(x, x_ind_row, x_ind_col)
      }
      stop(e)
    }
  )

  # add the loci table to the object
  attr(res, "loci") <- show_loci(x)
  class(res) <- c("gt_pcadapt", class(res))
  # Return the result
  return(res)
}

# Stop if some loci have the same genotype in every individual (monomorphic,
# or heterozygous in every individual, as heavy imputation can make them):
# pcadapt can not test them. It reads all the genotypes (one pass), so it is
# only called after snp_pcadapt() has failed.
# @param x a gen_tibble
# @param ind_row,ind_col the rows and columns of its backing matrix
# @keywords internal
# @noRd
stopifnot_variable_loci <- function(x, ind_row, ind_col) {
  stats <- bigstatsr::big_colstats(
    .gt_get_fbm(x),
    ind.row = ind_row,
    ind.col = ind_col
  )
  constant <- which(stats$var < 1e-12)
  if (length(constant) > 0) {
    stop(
      "gt_pcadapt() can not test loci with the same genotype in every ",
      "individual, and ", length(constant), " loci are monomorphic or ",
      "heterozygous in every individual (e.g. ",
      paste(utils::head(loci_names(x)[constant], 3), collapse = ", "),
      "), as heavy imputation can make them: remove them (e.g. with ",
      "select_loci()), or remove the loci with many missing genotypes ",
      "before imputing",
      call. = FALSE
    )
  }
}
