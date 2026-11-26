#' Compute the Genomic Relationship Matrix for a `gen_tibble` object
#'
#' This function computes the Genomic Relationship Matrix (GRM), with one of
#' two estimators.
#'
#' With `method = "WG17"` (the default), the GRM is estimated by computing the
#' pairwise kinship coefficients (coancestries) between all pairs of
#' individuals from a matrix of Allele Sharing following the approach of Weir
#' and Goudet 2017 based on beta estimators). The GRM is twice the coancestry
#' matrix (e.g. as estimated by `hierfstat::beta.dosage()` with `inb=FALSE`).
#' Every locus has the same weight, and relatedness is measured relative to
#' the average pair of individuals.
#'
#' With `method = "plink"`, the GRM is the variance-standardised relationship
#' matrix of PLINK 1.9 (`plink --make-rel`), the matrix whose eigenvectors are
#' the principal components of `plink --pca` and of [gt_pca()]:
#' \eqn{G_{ij} = \sum_l z_{il} z_{jl} / n_{ij}}, with standardised genotypes
#' \eqn{z_{il} = (g_{il} - 2p_l) / \sqrt{2p_l(1-p_l)}} (so that rare variants
#' count more), missing genotypes set to 0, and \eqn{n_{ij}} the number of
#' loci genotyped in both individuals. Pairs with no locus genotyped in both
#' are `NaN`. See [gt_pca()] for details.
#'
#' @param x a `gen_tibble` object.
#' @param allele_sharing_mat optional, the matrix of Allele Sharing returned by
#'   [pairwise_allele_sharing()] with `as_matrix=TRUE`. As a number of
#'   statistics can be derived from the Allele Sharing matrix, it is sometimes
#'   more efficient to pre-compute this matrix.
#' @param block_size the size of the blocks to use for the computation of the
#'   allele sharing matrix (with `method = "plink"`, blocks are chosen as in
#'   [gt_pca()]).
#' @param n_cores number of cores to be used for the computation of the
#'   allele sharing matrix, it defaults to [bigstatsr::nb_cores()]
#' @param method the estimator: `"WG17"` (Weir and Goudet 2017, the default)
#'   or `"plink"` (PLINK 1.9's variance-standardised matrix).
#' @returns a matrix of GR between all pairs of individuals
#' @export
#' @seealso [hierfstat::beta.dosage()], [gt_pca()]
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # Compute the GRM from the allele sharing matrix
#' example_gt %>% pairwise_grm()
#'
#' # To calculate using a precomputed allele sharing matrix, use:
#' allele_sharing <- example_gt %>% pairwise_allele_sharing(as_matrix = TRUE)
#' example_gt %>% pairwise_grm(allele_sharing_mat = allele_sharing)
#'
#' # The variance-standardised matrix of PLINK 1.9, as used by gt_pca()
#' example_gt %>% pairwise_grm(method = "plink")
pairwise_grm <- function(
  x,
  allele_sharing_mat = NULL,
  block_size = bigstatsr::block_size(nrow(x)),
  n_cores = bigstatsr::nb_cores(),
  method = c("WG17", "plink")
) {
  # only the autosomal loci (and pseudo-autosomal regions) are used
  # (see gt_chromosome_types())
  x <- .gt_autosomes_only(x)
  stopifnot_diploid(x)
  method <- match.arg(method)
  if (method == "plink") {
    if (!is.null(allele_sharing_mat)) {
      stop("'allele_sharing_mat' is only used with method = 'WG17'")
    }
    scaling <- pca_scaling(x, fun_scaling = NULL, n_cores = n_cores)
    grm <- grm_plink(x,
      center = scaling$center, scale = scaling$scale, n_cores = n_cores,
      block_size = grm_block_size(nrow(x)),
      use_blas = grm_use_blas("auto", n_cores)
    )
    # the count of missing genotypes is only used by gt_pca()
    attr(grm, "n_missing") <- NULL
    dimnames(grm) <- list(x$id, x$id)
    # add a class to indicate this is a pairwise matrix
    class(grm) <- c("pairwise_matrix", class(grm))
    return(grm)
  }
  if (is.null(allele_sharing_mat)) {
    allele_sharing_mat <- pairwise_allele_sharing(
      x,
      as_matrix = TRUE,
      block_size = block_size,
      n_cores = n_cores
    )
  }
  mii <- diag(allele_sharing_mat)
  diag(allele_sharing_mat) <- NA
  mb <- mean(allele_sharing_mat, na.rm = TRUE)
  diag(allele_sharing_mat) <- mii
  # estimate the pairwise kinship (coancestries)
  coa <- (allele_sharing_mat - mb) / (1 - mb)
  # add a class to indicate this is a pairwise matrix
  class(coa) <- c("pairwise_matrix", class(coa))
  return(coa * 2)
}
