#' Compute the Identity by State Matrix for a `gen_tibble` object
#'
#' This function computes the IBS matrix.
#'
#' Note that monomorphic sites are currently counted. Should we filter
#' them beforehand? What does plink do?
#' @param x a `gen_tibble` object.
#' @param as_counts whether the counts of similar alleles, rather than the proportion,
#' should be returned (FALSE by default).
#' @param block_size maximum number of loci read at once. More loci should improve speed,
#' but will tax memory.
#' @returns a list of two [bigstatsr::FBM] matrices, one of counts of IBS by alleles (i.e. 2*n loci),
#' and one of valid alleles (i.e. 2*n_loci - 2*missing_loci)
#' @export
gt_ibs <- function(x,
                      as_counts = FALSE,
                      block_size = bigstatsr::block_size(length(show_loci_names(x)))*2) {
  X <- attr(x$genotypes,"bigsnp") # convenient pointer
  x_ind_col <- show_loci(x)$big_index
  x_ind_row <- vctrs::vec_data(x$genotypes)
  ibs_matrix <- snp_ibs(X$genotypes,
            ind.row = x_ind_row,
            ind.col = x_ind_col,
            as.counts = as_counts,
            block.size = block_size)
  #rownames(ibs_matrix$counts_IBS) <- colnames(ibs_matrix$counts_IBS) <- x$id
  ibs_matrix
}
