#' Compute the KING-robust Matrix for a a `gen_tibble` object
#'
#' This function computes the KING-robust estimator of kinship.
#'
#' Note that monomorphic sites are currently considered. Remove monomorphic
#' sites before running pairwise_king if this is a concern.
#' @param x a `gen_tibble` object.
#' @param as_matrix boolean, determining whether the results should be a square
#'   symmetrical matrix (TRUE), or a tidied tibble (FALSE, the default)
#' @param block_size maximum number of loci read at once. More loci should
#'   improve speed, but will tax memory.
#' @returns a square symmetrical matrix of relationship coefficients between
#'   individuals if `as_matrix` is TRUE, or a tidied tibble of coefficients if
#'   `as_matrix` is FALSE.
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # Compute the KING-robust matrix
#' pairwise_king(example_gt, as_matrix = TRUE)
#'
#' # Or return a tidy tibble
#' pairwise_king(example_gt, as_matrix = FALSE)
#'
#' # Adjust block_size
#' pairwise_king(example_gt, block_size = 2)
#'
pairwise_king <- function(
    x,
    as_matrix = FALSE,
    block_size = bigstatsr::block_size(nrow(x))) {
  # nolint
  X <- attr(x$genotypes, "bigsnp") # convenient pointer #nolint
  x_ind_col <- show_loci(x)$big_index
  x_ind_row <- vctrs::vec_data(x$genotypes)
  king_matrix <- snp_king(
    X$genotypes,
    ind.row = x_ind_row,
    ind.col = x_ind_col,
    block.size = block_size
  )
  dimnames(king_matrix) <- list(x$id, x$id)
  if (as_matrix) {
    return(king_matrix)
  } else {
    return(tidy_dist_matrix(king_matrix))
  }
}
