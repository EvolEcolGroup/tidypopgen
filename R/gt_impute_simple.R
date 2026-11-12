#' Simple imputation based on allele frequencies
#'
#' This function provides a very simple imputation algorithm for `gen_tibble`
#' objects by using the mode, mean or sampling from the allele frequencies. Each
#' locus is imputed independently (and thus linkage information is ignored).
#'
#' This function is a wrapper around [bigsnpr::snp_fastImputeSimple()], which
#' imputes the whole file backing matrix (FBM). The allele frequencies used to
#' impute each locus are therefore estimated from all the individuals in the
#' FBM, so the `gen_tibble` must contain all of them (in any order). If some
#' individuals have been removed (e.g. with `filter()`), use
#' [gt_update_backingfile()] to resave the FBM with only the individuals in the
#' `gen_tibble` before imputing.
#'
#' As each locus is imputed independently, the `gen_tibble` does not need to
#' contain all the loci in the FBM: the loci that have been removed (e.g. with
#' [select_loci()]) are also imputed, but they do not affect the imputation of
#' the loci in the `gen_tibble`. If only a few of the loci in the FBM are kept,
#' resaving the FBM with [gt_update_backingfile()] first avoids the time spent
#' imputing loci that are not used.
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
#' @seealso [bigsnpr::snp_fastImputeSimple()] which this function wraps.
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
  n_cores = 1
) {
  method <- match.arg(method)

  if (n_cores > 1) {
    # Remove checking for two levels of parallelism
    .old_opt <- getOption("bigstatsr.check.parallel.blas", TRUE)
    options(bigstatsr.check.parallel.blas = FALSE)
    on.exit(options(bigstatsr.check.parallel.blas = .old_opt), add = TRUE)
  }

  # bigsnpr imputes all the rows and columns of the FBM. Loci are imputed
  # independently, so loci that are not in the gen_tibble are imputed but do
  # not change the loci that are; the allele frequencies, on the other hand,
  # are estimated from all the rows, so the gen_tibble must contain each
  # individual of the FBM exactly once (in any order)
  fbm_rows <- .gt_fbm_rows(x)
  if (
    length(fbm_rows) != nrow(attr(x$genotypes, "fbm")) ||
      anyDuplicated(fbm_rows) > 0
  ) {
    stop(
      "The individuals in the gen_tibble do not match the rows of the file ",
      "backing matrix. The imputation estimates the allele frequencies from ",
      "all the individuals in the file backing matrix, not only those in the ",
      "gen_tibble, so the file backing matrix has to be resaved first: use ",
      "gt_update_backingfile() to create a new file backing matrix with only ",
      "the individuals in the gen_tibble, and then impute it."
    )
  }

  if (gt_has_imputed(x)) {
    stop("object x is already imputed; use `gt_set_imputed(x, set = TRUE)`")
  }
  # imputed values would be stored in the bytes used by phased genotypes
  if (is_phased(x)) {
    stop(
      "phased genotypes can not be imputed; to impute the genotypes ",
      "without their phase, read the data with gen_tibble()"
    )
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

  attr(x$genotypes, "fbm") <- bigsnpr::snp_fastImputeSimple(
    attr(x$genotypes, "fbm"),
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
