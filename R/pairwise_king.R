#' Compute the KING-robust Matrix for a a `gen_tibble` object
#'
#' This function computes the KING-robust estimator of kinship.
#'
#' Note that monomorphic sites are currently considered. What does PLINK do???
#' @param x a `gen_tibble` object.
#' @param block_size maximum number of loci read at once. More loci should improve speed,
#' but will tax memory.
#' @export
gt_king <- function(x,
                      block_size = bigstatsr::block_size(length(show_loci_names(x)))) {

  X <- attr(x$genotypes,"bigsnp") # convenient pointer
  x_ind_col <- show_loci(x)$big_index
  x_ind_row <- vctrs::vec_data(x$genotypes)
  king_matrix <- snp_king(X$genotypes,
            ind.row = x_ind_row,
            ind.col = x_ind_col,
            block.size = block_size)
  rownames(king_matrix) <- colnames(king_matrix) <- x$id
  king_matrix
}
