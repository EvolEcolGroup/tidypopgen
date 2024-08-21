#' pcadapt analysis on a `gen_tibble` object
#'
#' pcadapt is an algorithm that detects genetic markers under selection. It is based on the
#' principal component analysis (PCA) of the genotypes of the individuals.
#' The method is described in [Luu et al. (2017)](https://doi.org/10.1534/genetics.116.195214),
#' See the R package `pcadapt`, which provides extensive
#' documentation and examples.
#'
#' Internally, this function uses the `snp_pcadapt` function from the `bigsnpr` package.
#' @param x A `gen_tibble` object.
#' @param pca a [`gt_pca`] object, as returned by `gt_pca_partialSVD()` or `gt_pca_randomSVD()`.
#' @param k Number of principal components to use in the analysis.
#' @param n_cores Number of cores to use.
#' @returns An object of subclass `gt_pcadapt`, a subclass of `mhtest`.
#' @export

gt_pcadapt <- function(x, pca, k, n_cores = 1) {
  stopifnot_gen_tibble(x)
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
  if (gt_has_imputed(x) && gt_uses_imputed(x)==FALSE){
    gt_set_imputed(x, set = TRUE)
    on.exit(gt_set_imputed(x, set = FALSE))
  }

  # Run the analysis
  res <- bigsnpr::snp_pcadapt(
    G = .gt_get_bigsnp(x)$genotypes,
    U.row = pca$u[, 1:k, drop = FALSE],
    ind.row = .gt_bigsnp_rows(x),
    ind.col = .gt_bigsnp_cols(x),
    ncores = n_cores
  )

  # add the loci table to the object
  attr(res, "loci") <- show_loci(x)
  class(res) <- c("gt_pcadapt", class(res))
  # Return the result
  return(res)
}
