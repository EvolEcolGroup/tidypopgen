#' Compute the KING-robust Matrix for a bigSNP object
#'
#' This function computes the KING-robust estimator of kinship. Missing genotypes will have to be
#' imputed before using this function. Computing the IBS of 1k individuals for
#' a million markers will take over 10 mins, so the function will not scale to
#' huge datasets. But for smaller dataset, it allows estimating the IBS directly
#' in R.
#' @param X a [bigstatsr::FBM.code256] matrix (as found in the `genotypes`
#' slot of a [bigsnpr::bigSNP] object).
#' @param row.ind An optional vector of the row indices that are used.
#' If not specified, all rows are used. Don't use negative indices.
#' @param col.ind An optional vector of the column indices that are used. If not
#'  specified, all columns are used. Don't use negative indices.
#' @param block.size maximum number of columns read at once.

#' @export
snp_king_r <- function(X,
                      row.ind = bigstatsr::rows_along(X),
                      col.ind = bigstatsr::cols_along(X),
                      block.size = bigstatsr::block_size(nrow(X))) {
  # check that the dataset has no NAs
  if (any(bigstatsr::big_counts(X)[4,]>0)){
    stop("NAs are not allowed, impute the missing genotypes")
  }
  # KING numerator for a matrix slice
  king_numerator <- function(X, ind){
    X_sub <- X[,ind] # get a slice of the matrix
    X_sub_t <- t(X_sub) # precompute the transpose of that matrix (as we use it several times)

    2*((X_sub==1) %*% t(X_sub==1) - 2*((X_sub==0) %*% t(X_sub==2) + (X_sub==2) %*% t(X_sub==0)) )
  }
  king_num_matrix <- bigstatsr::big_apply(X,
                                            a.FUN = king_numerator,
                                            ind = col.ind,
                                            a.combine = "plus",
                                            block.size = block.size)
  # now we compute the denominator
  # sum of how many genotypes are 1
  row_sums_1 <- bigstatsr::big_counts(X,byrow=TRUE)[2,]


  king_denominator = matrix(rep(row_sums_1, nrow(X)), nrow = nrow(X), byrow = T) +
    matrix(rep(row_sums_1, nrow(X)), nrow = nrow(X), byrow = F)

  king_num_matrix / king_denominator
}
