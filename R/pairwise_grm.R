#' Compute the Genomic Relationship Matrix for a `gen_tibble` object
#'
#' This function computes the Genomic Relationship Matrix (GRM). This is
#' estimated by computing the pairwise kinship coefficients (coancestries)
#' between all pairs of individuals from a matrix of Allele Sharing following
#' the approach of Weir and Goudet 2017 based on beta estimators).
#'
#' The GRM is twice the coancestry matrix (e.g. as estimated by
#' `hierfstat::beta.dosage()` with `inb=FALSE`).
#'
#' @param x a `gen_tibble` object.
#' @param allele_sharing_mat optional, the matrix of Allele Sharing returned by
#'   [pairwise_allele_sharing()] with `as_matrix=TRUE`. As a number of
#'   statistics can be derived from the Allele Sharing matrix, it it sometimes
#'   more efficient to pre-compute this matrix.
#' @param block_size the size of the blocks to use for the computation of the
#'   allele sharing matrix.
#' @returns a matrix of GR between all pairs of individuals
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # Compute the GRM from the allele sharing matrix
#' example_gt %>% pairwise_grm()
#'
#' # To calculate using a precomputed allele sharing matrix, use:
#' allele_sharing <- example_gt %>% pairwise_allele_sharing(as_matrix = TRUE)
#' example_gt %>% pairwise_grm(allele_sharing_mat = allele_sharing)
pairwise_grm <- function(
    x,
    allele_sharing_mat = NULL,
    block_size = bigstatsr::block_size(nrow(x))) {
  if (is.null(allele_sharing_mat)) {
    allele_sharing_mat <- pairwise_allele_sharing(
      x,
      as_matrix = TRUE,
      block_size = block_size
    )
  }
  mii <- diag(allele_sharing_mat)
  diag(allele_sharing_mat) <- NA
  mb <- mean(allele_sharing_mat, na.rm = TRUE)
  diag(allele_sharing_mat) <- mii
  # estimate the pairwise kinship (coancestries)
  coa <- (allele_sharing_mat - mb) / (1 - mb)
  return(coa * 2)
}
