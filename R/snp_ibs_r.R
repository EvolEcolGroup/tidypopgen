#' Compute the Identity by State Matrix for a bigSNP object
#'
#' This function computes the IBS matrix. Missing genotypes will have to be
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
#' @param as.counts whether the counts of similar alleles, rather than the proportion,
#' should be returned (FALSE by default).
#' @param block.size maximum number of columns read at once.
#' @export
snp_ibs_r <- function(X,
                      row.ind = bigstatsr::rows_along(X),
                      col.ind = bigstatsr::cols_along(X),
                      as.counts = FALSE,
                      block.size = bigstatsr::block_size(nrow(X))) {
  # check that the dataset has no NAs
  if (any(bigstatsr::big_counts(X)[4,]>0)){
    stop("NAs are not allowed, impute the missing genotypes")
  }
  # IBS (as counts) for a matrix slice
  ibs_counts <- function(X, ind){
    X_sub <- X[,ind] # get a slice of the matrix
    X_sub_t <- t(X_sub)
    2 * ( (X_sub==2) %*% (X_sub_t==2) + (X_sub==1) %*% (X_sub_t==1) + (X_sub==0) %*% (X_sub_t==0) ) + #prop. of loci with same genotypes
      ( (X_sub==1) %*% (X_sub_t==0 | X_sub_t==2) + (X_sub==0 | X_sub==2) %*% (X_sub_t==1) ) #prop. of loci with one allele shared
  }
  ibs_counts_matrix <- bigstatsr::big_apply(X,
                                            a.FUN = ibs_counts,
                                            ind = col.ind,
                                            a.combine = "plus",
                                            block.size = block.size)
  if (as.counts){
    return(ibs_counts_matrix)
  } else {
    ibs_counts_matrix/(2*length(col.ind))
  }
}
