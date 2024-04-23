#' Compute the Identity by State Matrix for a `gen_tibble` object
#'
#' This function computes the IBS matrix.
#'
#' Note that monomorphic sites are currently counted. Should we filter
#' them beforehand? What does plink do?
#' @param x a `gen_tibble` object.
#' @param type one of "proportion" (equivalent to "ibs" in PLINK), "adjusted_counts" ("distance" in PLINK),
#' and "raw_counts" (the counts of identical alleles and non-missing alleles, from which the two other quantities are computed)
#' @param block_size maximum number of loci read at once. More loci should improve speed,
#' but will tax memory.
#' @returns a [bigstatsr::FBM] of proportion or adjusted counts, or a list of
#' two [bigstatsr::FBM] matrices, one of counts of IBS by alleles,
#' and one of number of valid alleles (i.e. 2*n_loci - 2*missing_loci)
#' @export
gt_ibs <- function(x,
                      type = c("proportion","adjusted_counts","raw_counts"),
                      block_size = bigstatsr::block_size(count_loci(x))) {
  type <- match.arg(type)
  X <- attr(x$genotypes,"bigsnp") # convenient pointer
  x_ind_col <- show_loci(x)$big_index
  x_ind_row <- vctrs::vec_data(x$genotypes)
  ibs_matrix <- snp_ibs(X$genotypes,
            ind.row = x_ind_row,
            ind.col = x_ind_col,
            type = type,
            block.size = block_size)
  if (inherits(ibs_matrix,"matrix")){
    dimnames(ibs_matrix)<-list(x$id, x$id)
  } else { # else if we have a list of two count matrices
    attr(ibs_matrix[[1]],"indiv_names") <- x$id
    attr(ibs_matrix[[2]],"indiv_names") <- x$id
  }
  ibs_matrix
}
