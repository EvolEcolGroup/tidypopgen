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
#' The genotypes tested are those that the PCA was computed on:
#' - with a PCA of genotypes with missing values (the default method of
#'   [gt_pca()], which compares individuals over the loci genotyped in both,
#'   as PLINK does), the genotypes of `x` with their missing values: each
#'   locus is tested on the individuals genotyped at it, as in the `pcadapt`
#'   package. Imputed genotypes are not used.
#' - with a PCA of complete or imputed genotypes, the genotypes of `x` must be
#'   complete or imputed (see [gt_impute_simple()]; imputed genotypes are
#'   used automatically). Missing genotypes give an error; as in `bigsnpr`,
#'   they are looked for in a sample of the genotypes (all the individuals,
#'   at evenly spaced loci), and any that go unnoticed are left out of the
#'   test of their locus.
#'
#' Loci with the same genotype in every individual genotyped, or genotyped
#' in fewer than 3 individuals, can not be tested, and give an error.
#' @references Luu, K., Bazin, E., Blum, M. G. B., & François, O. (2017).
#'   pcadapt: an R package for genome scans for selection based on principal
#'   component analysis. Molecular Ecology Resources, 17(1), 67–77.
#' @param x A `gen_tibble` object.
#' @param pca a [`gt_pca`] object, as returned by [gt_pca()], of the
#'   individuals of `x` (see Details for missing genotypes).
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
#' # Remove monomorphic loci
#' lobsters <- lobsters %>% select_loci_if(loci_maf(genotypes) > 0)
#'
#' # Create PCA object (the lobsters have some missing genotypes, which do not
#' # need to be imputed)
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
  # the genotypes to test: those the PCA was computed on
  pca_missing <- isTRUE(attr(pca, "missing"))
  if (pca_missing) {
    # a PCA of genotypes with missing values: test them (without imputed
    # genotypes); each locus is then tested on the individuals genotyped at
    # it, as in the pcadapt package
    if (gt_has_imputed(x) && gt_uses_imputed(x)) {
      gt_set_imputed(x, set = FALSE)
      on.exit(gt_set_imputed(x, set = TRUE), add = TRUE)
    }
  } else if (gt_has_imputed(x) && gt_uses_imputed(x) == FALSE) {
    # a PCA of complete or imputed genotypes: use the imputed genotypes
    gt_set_imputed(x, set = TRUE)
    on.exit(gt_set_imputed(x, set = FALSE), add = TRUE)
  }
  # bigsnpr's check for missing values samples the whole backing matrix
  # (including individuals and loci that x does not use), so turn it off
  # (keeping bigstatsr's check of n_cores); with a PCA of complete or imputed
  # genotypes, check a sample of the genotypes of x instead
  local_no_bigstatsr_checks(n_cores = n_cores)
  if (!pca_missing) {
    stopifnot_no_missing(x, "gt_pcadapt()", n_cores = n_cores)
  }

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
      # can't have missing values in 'U'" (#179), or, with k = 1, as an
      # error of its genomic control ("invalid 'xmin' value"): look for them
      # (only after an error) and say which they are
      stopifnot_variable_loci(x, x_ind_row, x_ind_col)
      stop(e)
    }
  )

  # add the loci table to the object
  attr(res, "loci") <- show_loci(x)
  class(res) <- c("gt_pcadapt", class(res))
  # Return the result
  return(res)
}

# Stop if some loci have the same genotype in every individual genotyped
# (monomorphic, or heterozygous in every individual, as heavy imputation can
# make them), or are genotyped in fewer than 3 individuals: pcadapt can not
# test them (their t-scores are NA in bigsnpr). It reads all the genotypes
# (one pass), so it is only called after snp_pcadapt() has failed.
# @param x a gen_tibble
# @param ind_row,ind_col the rows and columns of its backing matrix
# @keywords internal
# @noRd
stopifnot_variable_loci <- function(x, ind_row, ind_col) {
  # the counts of each genotype (0, 1 and 2, and missing, NA) at each locus
  counts <- bigstatsr::big_counts(
    .gt_get_fbm(x),
    ind.row = ind_row,
    ind.col = ind_col
  )
  counts <- counts[!is.na(rownames(counts)), , drop = FALSE]
  # loci with a single genotype among the individuals genotyped there, or
  # with fewer than 3 of them (a perfect fit with 2)
  constant <- which(colSums(counts > 0) < 2 | colSums(counts) < 3)
  if (length(constant) > 0) {
    stop(
      "gt_pcadapt() can not test loci with the same genotype in every ",
      "individual, and ", length(constant), " loci are monomorphic or ",
      "heterozygous in every individual (genotyped), or genotyped in fewer ",
      "than 3 individuals (e.g. ",
      paste(utils::head(loci_names(x)[constant], 3), collapse = ", "),
      "): remove them (e.g. with select_loci()); heavy imputation can make ",
      "loci heterozygous in every individual, so remove the loci with many ",
      "missing genotypes before imputing",
      call. = FALSE
    )
  }
}
