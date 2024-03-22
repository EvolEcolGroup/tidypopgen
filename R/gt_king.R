#' Compute the KING-robust Matrix for a a `gen_tibble` object
#'
#' This function computes the KING-robust estimator of kinship. Missing genotypes will have to be
#' imputed before using this function. Computing the IBS of 1k individuals for
#' a million markers will take over 10 mins, so the function will not scale to
#' huge datasets. But for smaller dataset, it allows estimating the IBS directly
#' in R.
#'
#' Note that monomorphic sites are currently considered. It would be better to filter
#' them beforehand.
#' @param x a `gen_tibble` object.
#' @param as_counts whether the counts of similar alleles, rather than the proportion,
#' should be returned (FALSE by default).
#' @param block_size maximum number of loci read at once. More loci should improve speed,
#' but will tax memory.
#' @export
gt_king <- function(x,
                      as_counts = FALSE,
                      block_size = bigstatsr::block_size(length(show_loci_names(x)))) {

  X <- attr(x$genotypes,"bigsnp") # convenient pointer
  x_ind_col <- show_loci(x)$big_index
  x_ind_row <- vctrs::vec_data(x$genotypes)
  king_matrix <- snp_king_r(X$genotypes,
            row.ind = x_ind_row,
            col.ind = x_ind_col,
            block.size = block_size)
  rownames(king_matrix) <- colnames(king_matrix) <- x$id
  king_matrix
}
