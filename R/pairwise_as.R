#' Compute the Pairwise Allele Sharing Matrix for a `gen_tibble` object
#'
#' This function computes the Allele Sharing matrix.
#' Estimates Allele Sharing (matching in `hierfstat`)) between pairs of individuals
#' (for each locus, gives 1 if the two individuals are homozygous
#' for the same allele, 0 if they are homozygous for a different allele, and 1/2 if at least one individual
#' is heterozygous. Matching is the average of these 0, 1/2 and 1s)
#'
#' @param x a `gen_tibble` object.
#' @param block_size maximum number of loci read at once. More loci should improve speed,
#' but will tax memory.
#' @returns a matrix of allele sharing between all pairs of individuals
#' @export
pairwise_as <- function(x,
                      block_size = bigstatsr::block_size(count_loci(x))) {
  X <- attr(x$genotypes,"bigsnp") # convenient pointer
  x_ind_col <- .gt_bigsnp_cols(x)
  x_ind_row <- .gt_bigsnp_rows(x)
  as_matrix <- snp_as(X$genotypes,
            ind.row = x_ind_row,
            ind.col = x_ind_col,
            block.size = block_size)
  dimnames(as_matrix)<-list(x$id, x$id)
  as_matrix
}
