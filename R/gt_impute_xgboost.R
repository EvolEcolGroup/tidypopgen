#' Imputation based XGBoost
#'
#' This function provides a simple imputation algorithm for `gen_tibble`
#' objects based on local XGBoost models.
#'
#' This function is a wrapper around [bigsnpr::snp_fastImpute()].
#'
#' @param x a [gen_tibble] with missing data
#' @param alpha Type-I error for testing correlations. Default is `1e-4`.
#' @param size Number of neighbor SNPs to be possibly included in the model
#'   imputing this particular SNP. Default is `200`.
#' @param p_train Proportion of non missing genotypes that are used for training
#'   the imputation model while the rest is used to assess the accuracy of
#'   this imputation model. Default is `0.8`.
#' @param n_cor Number of rows that are used to estimate correlations.
#'   Default uses them all.
#' @param seed An integer, for reproducibility. Default doesn't use seeds.
#' @param n_cores the number of cores to be used
#' @returns a [gen_tibble] with imputed genotypes
#' @export

gt_impute_xgboost <- function(
    x,
    alpha = 1e-4,
    size = 200,
    p_train = 0.8,
    n_cor = nrow(x),
    seed = NA,
    n_cores = 1) {
  method <- match.arg(method)

  if (n_cores > 1) {
    # Remove checking for two levels of parallelism
    options(bigstatsr.check.parallel.blas = FALSE)
    on.exit(options(bigstatsr.check.parallel.blas = TRUE))
  }

  if (nrow(x) != nrow(attr(x$genotypes, "bigsnp")$genotypes)) {
    stop(
      "The number of individuals in the gen_tibble does not match the",
      " number of rows in the file backing matrix. Before imputing, use",
      " gt_update_backingfile to update your file backing matrix."
    )
  }

  if (
    !all.equal(attr(x$genotypes, "bigsnp")$genotypes$code256, bigsnpr::CODE_012)
  ) {
    # nolint start
    if (
      all.equal(
        attr(x$genotypes, "bigsnp")$genotypes$code256,
        bigsnpr::CODE_IMPUTE_PRED
      ) ||
        all.equal(
          attr(x$genotypes, "bigsnp")$genotypes$code256,
          bigsnpr::CODE_DOSAGE
        )
    ) {
      # nolint end
      stop("object x is already imputed")
    } else {
      stop("object x uses a code256 that is not compatible with imputation")
    }
  }

  attr(x$genotypes, "bigsnp")$genotypes <- bigsnpr::snp_fastImputeSimple(
    attr(x$genotypes, "bigsnp")$genotypes, # this needs subsetting
    infos.chr = show_loci(x)$chr_int, # check this is correct
                               alpha = alpha,
                               size = size,
                               p.train = p_train,
                               n.cor = n_cor,
                               seed = seed,
                               ncores = n_cores)

  attr(x$genotypes, "imputed") <- "simple"
  gt_set_imputed(x, set = FALSE)
  x
}
