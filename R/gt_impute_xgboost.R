#' Imputation based XGBoost
#'
#' This function provides a simple imputation algorithm for `gen_tibble`
#' objects based on local XGBoost models.
#'
#' This function is a wrapper around [bigsnpr::snp_fastImpute()]. The error
#' rates from the xgboost, if appended, can be retrieved with
#' `attr(x$genotypes, "imputed_errors")` where `x` is the `gen_tibble`.
#'
#' @param x a [gen_tibble] with missing data
#' @param alpha Type-I error for testing correlations. Default is `1e-4`.
#' @param size Number of neighbour SNPs to be possibly included in the model
#'   imputing this particular SNP. Default is `200`.
#' @param p_train Proportion of non missing genotypes that are used for training
#'   the imputation model while the rest is used to assess the accuracy of
#'   this imputation model. Default is `0.8`.
#' @param n_cor Number of rows that are used to estimate correlations.
#'   Default uses them all.
#' @param seed An integer, for reproducibility. Default doesn't use seeds.
#' @param n_cores the number of cores to be used
#' @param append_error boolean, should the xgboost error estimates be appended
#' as an attribute to the genotype column of
#' the gen_tibble. If TRUE (the default), a matrix of two rows (the number
#' of missing values, and the error estimate) and as many columns as the number
#' of loci will be appended to the gen_tibble.
#' attr(missing_gt$genotypes, "imputed_errors")
#' @returns a [gen_tibble] with imputed genotypes
#' @export
#' @seealso [bigsnpr::snp_fastImpute()] which this function wraps.
#' @examplesIf rlang::is_installed(c("xgboost", "RhpcBLASctl", "data.table"))
#' \dontshow{
#' data.table::setDTthreads(2)
#' RhpcBLASctl::blas_set_num_threads(2)
#' RhpcBLASctl::omp_set_num_threads(2)
#' }
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # Impute the gen_tibble
#' example_gt <- example_gt %>% gt_impute_xgboost()
#'
#' # And we can check it has been imputed
#' example_gt %>% gt_has_imputed()
gt_impute_xgboost <- function(
    x,
    alpha = 1e-4,
    size = 200,
    p_train = 0.8,
    n_cor = nrow(x),
    seed = NA,
    n_cores = 1,
    append_error = TRUE) {
  if (n_cores > 1) {
    # Remove checking for two levels of parallelism
    .old_opt <- getOption("bigstatsr.check.parallel.blas", TRUE)
    options(bigstatsr.check.parallel.blas = FALSE)
    on.exit(options(bigstatsr.check.parallel.blas = .old_opt), add = TRUE)
  }

  if (nrow(x) != nrow(attr(x$genotypes, "fbm"))) {
    stop(
      "The number of individuals in the gen_tibble does not match the",
      " number of rows in the file backing matrix. Before imputing, use",
      " gt_update_backingfile to update your file backing matrix."
    )
  }

  if (count_loci(x) != ncol(attr(x$genotypes, "fbm"))) {
    stop(
      "The number of loci in the gen_tibble does not match the number of ",
      "columns in the file backing matrix. Before imputing, use ",
      "gt_update_backingfile to update your file backing matrix."
    )
  }

  if (gt_has_imputed(x)) {
    stop("object x is already imputed; use `gt_set_imputed(x, set = TRUE)`")
  }

  if (
    !identical(attr(x$genotypes, "fbm")$code256, bigsnpr::CODE_012)
  ) {
    # nolint start
    if (
      identical(
        attr(x$genotypes, "fbm")$code256,
        bigsnpr::CODE_IMPUTE_PRED
      ) ||
        identical(
          attr(x$genotypes, "fbm")$code256,
          bigsnpr::CODE_DOSAGE
        )
    ) {
      # nolint end
      stop("object x is already imputed, but attr(x, 'imputed') is null")
    } else {
      stop("object x uses a code256 that is not compatible with imputation")
    }
  }

  infos <- bigsnpr::snp_fastImpute(
    attr(x$genotypes, "fbm"), # this needs subsetting
    infos.chr = cast_chromosome_to_int(show_loci(x)$chromosome),
    alpha = alpha,
    size = size,
    p.train = p_train,
    n.cor = n_cor,
    seed = seed,
    ncores = n_cores
  )

  attr(x$genotypes, "imputed") <- "xgboost"
  if (append_error) {
    attr(x$genotypes, "imputed_errors") <- infos[]
  }
  # TODO add a message if append_error = FALSE - something like
  # “an FBM has been generated in the temp directory,
  # and this can be found here : ...”
  gt_set_imputed(x, set = FALSE)
  x
}
