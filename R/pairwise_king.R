#' Compute the KING-robust Matrix for a a `gen_tibble` object
#'
#' This function computes the KING-robust estimator of kinship.
#'
#' Note that monomorphic sites are currently considered. What does PLINK do???
#' @param x a `gen_tibble` object.
#' @param as_matrix boolean, determining whether the results should be a square symmetrical matrix (TRUE),
#' or a tidied tibble (FALSE, the default)
#' @param block_size maximum number of loci read at once. More loci should improve speed,
#' but will tax memory.
#' @export
gt_king <- function(x, as_matrix = FALSE,
                      block_size = bigstatsr::block_size(length(show_loci_names(x)))) {

  X <- attr(x$genotypes,"bigsnp") # convenient pointer
  x_ind_col <- show_loci(x)$big_index
  x_ind_row <- vctrs::vec_data(x$genotypes)
  king_matrix <- snp_king(X$genotypes,
            ind.row = x_ind_row,
            ind.col = x_ind_col,
            block.size = block_size)
  dimnames(king_matrix)<-list(x$id, x$id)
  if (as_matrix){
    return(king_matrix)
  } else {
    return(tidy_dist_matrix(king_matrix))
  }

}
