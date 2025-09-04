#' Simple imputation based on allele frequencies
#'
#' This function provides a very simple imputation algorithm for `gen_tibble`
#' objects by using the mode, mean or sampling from the allele frequencies. Each
#' locus is imputed independently (and thus linkage information is ignored).
#'
#' This function is a wrapper around [bigsnpr::snp_fastImputeSimple()].
#'
#' @param x a [gen_tibble] with missing data
#' @param method one of
#' - 'mode': the most frequent genotype
#' - 'mean0': the mean rounded to the nearest integer
#' - 'random': randomly sample a genotype based on the observed
#'    allele frequencies
#' @param n_cores the number of cores to be used
#' @returns a [gen_tibble] with imputed genotypes
#' @export
#' @examplesIf all(rlang::is_installed(c("RhpcBLASctl", "data.table")))
#' \dontshow{
#' data.table::setDTthreads(2)
#' RhpcBLASctl::blas_set_num_threads(2)
#' RhpcBLASctl::omp_set_num_threads(2)
#' }
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # Impute the gen_tibble
#' example_gt <- example_gt %>% gt_impute_simple()
#'
#' # And we can check it has been imputed
#' example_gt %>% gt_has_imputed()
gt_impute_simple <- function(
    x,
    method = c("mode", "mean0", "random"),
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

  if (gt_has_imputed(x)) {
    stop("object x is already imputed, use `gt_set_imputed(x, TRUE)`")
  }

  if (
    !identical(attr(x$genotypes, "bigsnp")$genotypes$code256, bigsnpr::CODE_012)
  ) {
    # nolint start
    if (
      identical(
        attr(x$genotypes, "bigsnp")$genotypes$code256,
        bigsnpr::CODE_IMPUTE_PRED
      ) ||
        identical(
          attr(x$genotypes, "bigsnp")$genotypes$code256,
          bigsnpr::CODE_DOSAGE
        )
    ) {
      # nolint end
      stop("object x is already imputed, but attr(x, 'imputed') is null")
    } else {
      stop("object x uses a code256 that is not compatible with imputation")
    }
  }

  attr(x$genotypes, "bigsnp")$genotypes <- bigsnpr::snp_fastImputeSimple(
    attr(x$genotypes, "bigsnp")$genotypes,
    method = method,
    ncores = n_cores
  )

  attr(x$genotypes, "imputed") <- "simple"
  # prioritise "gen_tbl" class over "sf"
  obj_class <- class(x)
  if ("sf" %in% obj_class) {
    obj_class <-
      c("gen_tbl", "sf", obj_class[!obj_class %in% c("gen_tbl", "sf")])
    class(x) <- obj_class
  }
  gt_set_imputed(x, set = FALSE)
  x
}
