#' Simple imputation based on allele frequencies
#'
#' This function provides a very simple imputation algorithm for `gen_tibble`
#' objects by using the mode, mean or sampling from the allele frequencies.
#' Each locus is imputed independently (and thus linkage information is ignored).
#'
#' This function is a wrapper around [bigsnpr::snp_fastImputeSimple()].
#'
#' @param x a [gen_tibble] with missing data
#' @param method one of
#' - 'mode': the most frequent genotype
#' - 'mean0': the mean rounded to the nearest integer
#' - 'random': randomly sample a genotype based on the observed allele frequencies
#' @param n_cores the number of cores to be used
#' @returns a [gen_tibble] with imputed genotypes
#' @export

gt_impute_simple <-function(x,
                            method = c("mode", "mean0", "random"),
                            n_cores = 1) {
  method <- match.arg(method)

  if (!all.equal(attr(x$genotypes,"bigsnp")$genotypes$code256,bigsnpr::CODE_012)){
    if (all.equal(attr(x$genotypes,"bigsnp")$genotypes$code256, bigsnpr::CODE_IMPUTE_PRED) |
        all.equal(attr(x$genotypes,"bigsnp")$genotypes$code256, bigsnpr::CODE_DOSAGE)){
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
  attr(x$genotypes,"imputed")<-"simple"
  gt_set_imputed(x, set= FALSE)
  x
}
