#' Compute the Pairwise Allele Sharing Matrix for a `gen_tibble` object
#'
#' This function computes the Allele Sharing matrix. Estimates Allele Sharing
#' (equivalent to the quantity estimated by `hierfstat::matching()`) between
#' pairs of individuals (for each locus, gives 1 if the two individuals are
#' homozygous for the same allele, 0 if they are homozygous for a different
#' allele, and 1/2 if at least one individual is heterozygous. Matching is the
#' average of these 0, 1/2 and 1s)
#'
#' @param x a `gen_tibble` object.
#' @param as_matrix boolean, determining whether the results should be a square
#'   symmetrical matrix (TRUE), or a tidied tibble (FALSE, the default)
#' @param block_size maximum number of loci read at once. More loci should
#'   improve speed, but will tax memory.
#' @returns a matrix of allele sharing between all pairs of individuals
#' @export
#' @seealso [hierfstat::matching()]
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # Compute allele sharing between individuals
#' example_gt %>% pairwise_allele_sharing(as_matrix = FALSE)
#'
#' # Alternatively, return as a tibble
#' example_gt %>% pairwise_allele_sharing(as_matrix = TRUE)
pairwise_allele_sharing <- function(
    x,
    as_matrix = FALSE,
    block_size = bigstatsr::block_size(nrow(x))) {
  # nolint
  X <- attr(x$genotypes, "fbm") # convenient pointer #nolint
  x_ind_col <- .gt_fbm_cols(x)
  x_ind_row <- .gt_fbm_rows(x)
  ashare_matrix <- snp_allele_sharing(
    X,
    ind.row = x_ind_row,
    ind.col = x_ind_col,
    block.size = block_size
  )
  dimnames(ashare_matrix) <- list(x$id, x$id)
  if (as_matrix) {
    return(ashare_matrix)
  } else {
    return(tidy_dist_matrix(ashare_matrix))
  }
}
