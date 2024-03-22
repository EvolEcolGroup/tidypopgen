#' Compute the KING-robust Matrix for a bigSNP object
#'
#' This function computes the KING-robust estimator of kinship.
#'
#' The results should be equivalent to PLINK, but that SHOULD be tested!!!!
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
                      block.size = bigstatsr::block_size(nrow(X))*3) {
  #check_args()

  n <- length(ind.row)
  # FBM matrix for the KING numerator
  K_numerator <- bigstatsr::FBM(n, n, init = 0)
  m <- length(ind.col)

  intervals <- CutBySize(m, block.size)

  # Preassign memory where we will store the slices of genotypes.
  # For efficiency, when we read in the slice,
  # we will immediately recode it into these3 matrices,
  # one per genotype, equivalent to X==0, X==1, and X==2 respectively
  X_part_temp0 <- matrix(0, n, max(intervals[, "size"]))
  X_part_temp1 <- matrix(0, n, max(intervals[, "size"]))
  X_part_temp2 <- matrix(0, n, max(intervals[, "size"]))

  for (j in bigstatsr::rows_along(intervals)) {
    ind <- seq2(intervals[j, ]) # this iteration indices
    ind.col.ind <- ind.col[ind] #subset given indices by the iteration indeces
    increment_king_numerator(K_numerator,
                         X_part_temp0,
                         X_part_temp1,
                         X_part_temp2,
                         X,
                         ind.row,
                         ind.col.ind)
  }
  # now we compute the denominator
  # sum of how many genotypes are 1
  row_sums_1 <- bigstatsr::big_counts(X,byrow=TRUE)[2,]

  # divide KING num by den for a set of rows 'ind'
  divide_king_sub <- function (K_num_sub, ind, row_sums_1){
    X_sub <- K_num_sub[,ind]
    n_rows <- nrow(K_num_sub)
    # denominator of KING for these columns
    K_den_sub <- K_den <- matrix(rep(row_sums_1[ind],nrow(K_num_sub)),nrow=nrow(K_num_sub), byrow=T)+
      matrix(rep(row_sums_1,length(ind)),nrow=nrow(K_num_sub), byrow=F)
    # the ratio
    X_sub/K_den_sub
  }
  # This works, but could be done with a FBM in c++ for added speed
  ibs_counts_matrix <- bigstatsr::big_apply(K_numerator,
                                            a.FUN = divide_king_sub,
                                            ind = seq_len(ncol(K_numerator)),
                                            row_sums_1 = row_sums_1,
                                            a.combine = "cbind",
                                            # in theory block size coudl be our previouse block size * 3
                                            block.size = bigstatsr::block_size(nrow(K_numerator)))
}
