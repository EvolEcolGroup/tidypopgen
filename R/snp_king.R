#' Compute the KING-robust Matrix for a bigSNP object
#'
#' This function computes the KING-robust estimator of kinship.
#'
#' The last step is not optimised yet, as it does the division of the num by the
#' den all in memory (on my TODO list...).
#'
#' @param X a [bigstatsr::FBM.code256] matrix (as found in the `genotypes`
#' slot of a [bigsnpr::bigSNP] object).
#' @param ind.row An optional vector of the row indices that are used.
#' If not specified, all rows are used. Don't use negative indices.
#' @param ind.col An optional vector of the column indices that are used. If not
#'  specified, all columns are used. Don't use negative indices.
#' @param block.size maximum number of columns read at once.

#' @export
snp_king <- function(X,
                      ind.row = bigstatsr::rows_along(X),
                      ind.col = bigstatsr::cols_along(X),
                      block.size = bigstatsr::block_size(nrow(X))*4) {
  if (!inherits(X,"FBM.code256")){
    stop ("X should be a FBM.code256 matrix")
  }

  n <- length(ind.row)
  # FBM matrix for the KING numerator
  K_numerator <- bigstatsr::FBM(n, n, init = 0)
  # FBM matrix of number of heterozygous loci for individual i (limited to sites not missing in j)
  N_Aa_i <- bigstatsr::FBM(n, n, init = 0)
  m <- length(ind.col)

  intervals <- CutBySize(m, block.size)

  # Preassign memory where we will store the slices of genotypes.
  # For efficiency, when we read in the slice,
  # we will immediately recode it into these3 matrices,
  # one per genotype, equivalent to X==0, X==1, and X==2 respectively
  X_0_part <- matrix(0, n, max(intervals[, "size"]))
  X_1_part <- matrix(0, n, max(intervals[, "size"]))
  X_2_part <- matrix(0, n, max(intervals[, "size"]))
  # valid loci (i.e. wiht genotype 0,1 or 2)
  X_valid_part <- matrix(0, n, max(intervals[, "size"]))

  for (j in bigstatsr::rows_along(intervals)) {
    ind <- seq2(intervals[j, ]) # this iteration indices
    ind.col.ind <- ind.col[ind] #subset given indices by the iteration indeces
    increment_king_numerator(K_numerator,
                             N_Aa_i,
                         X_0_part,
                         X_1_part,
                         X_2_part,
                         X_valid_part,
                         X,
                         ind.row,
                         ind.col.ind)
  }
  # transpose the counts of heterozygous sites to confirm that this is valid
  N_Aa_j <- bigstatsr::big_transpose(N_Aa_i)

  # divide KING num by den for a set of col 'ind'
  divide_king_sub <- function (K, ind, N_Aa_i, N_Aa_j){
    K_sub <- K[,ind]
    N_Aa_i_sub <- N_Aa_i[,ind]
    N_Aa_j_sub <- N_Aa_j[,ind]
    # denominator of KING for these columns
    K_sub/(2* pmin(N_Aa_i_sub,N_Aa_j_sub))+0.5-
      0.25*(N_Aa_i_sub+N_Aa_j_sub)/pmin(N_Aa_i_sub,N_Aa_j_sub)
  }
  # This works, but could be done with a FBM in c++ for added speed
  ibs_counts_matrix <- bigstatsr::big_apply(K_numerator,
                                            a.FUN = divide_king_sub,
                                            ind = seq_len(ncol(K_numerator)),
                                            N_Aa_i = N_Aa_i,
                                            N_Aa_j = N_Aa_j,
                                            a.combine = "cbind",
                                            block.size = bigstatsr::block_size(nrow(K_numerator)))
}
