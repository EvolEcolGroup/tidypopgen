#' Simple imputation based on allele frequencies
#'
#' This function provides a very simple imputation algorithm for `gen_tibble`
#' objects by using the mode, mean or sampling from the allele frequencies.
#' Each locus is imputed independently (and thus linkage information is ignored).
#' It is a wrapper around [bigsnpr::snp_fastImputeSimple()].
#'
#' @param x a [gen_tibble] with missing data
#' @param method one of
#' - 'median': the most frequent genotype
#' - 'mean0': the mean rounded to the nearest integer
#' - 'mean2': the mean rounded to 2 decimal places
#' - 'random': randomly sample a genotype based on the observed allele frequencies
#' @param n_cores the number of cores to be used
#' @returns a [gen_tibble] with imputed genotypes
#' @export

gt_impute_simple <-function(x,
                            method = c("mode", "mean0", "mean2", "random"),
                            n_cores = 1) {
  if (!identical(attr(x$genotypes,"bigsnp")$genotypes$code256,bigsnpr::CODE_012)){
    if (identical(attr(x$genotypes,"bigsnp")$genotypes$code256, bigsnpr::CODE_IMPUTE_PRED) |
        identical(attr(x$genotypes,"bigsnp")$genotypes$code256, bigsnpr::CODE_DOSAGE)){
      stop("object x is already imputed")
    } else {
      stop ("object x uses a code256 that is not compatible with imputation")
    }
  }

  attr(x$genotypes,"bigsnp")$genotypes <- bigsnpr::snp_fastImputeSimple(
    attr(x$genotypes,"bigsnp")$genotypes,
    method = method,
    ncores = n_cores
  )
  x
}