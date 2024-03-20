#' Compute the Identity by State Matrix for a `gen_tibble` object
#'
#' This function computes the IBS matrix. Missing genotypes will have to be
#' imputed before using this function. Computing the IBS of 1k individuals for
#' a million markers will take over 10 mins, so the function will not scale to
#' huge datasets. But for smaller dataset, it allows estimating the IBS directly
#' in R.
#'
#' Note that monomorphic sites are currently counted. It would be better to filter
#' them beforehand.
#' @param x a `gen_tibble` object.
#' @param as_counts whether the counts of similar alleles, rather than the proportion,
#' should be returned (FALSE by default).
#' @param block_size maximum number of loci read at once. More loci should improve speed,
#' but will tax memory.
#' @export
gt_ibs <- function(x,
                      as_counts = FALSE,
                      block_size = bigstatsr::block_size(length(show_loci_names(x)))) {

  X <- attr(x$genotypes,"bigsnp") # convenient pointer
  x_ind_col <- show_loci(x)$big_index
  x_ind_row <- vctrs::vec_data(x$genotypes)
  ibs_matrix <- snp_ibs_r(X$genotypes,
            row.ind = x_ind_row,
            col.ind = x_ind_col,
            as.counts = as_counts,
            block.size = block_size)
  rownames(ibs_matrix) <- colnames(ibs_matrix) <- x$id
  ibs_matrix
}
