#' Compute the Identity by State Matrix for a bigSNP object
#'
#' This function computes the IBS matrix.
#'
#' Note that monomorphic sites are currently counted. Should we filter
#' them beforehand? What does plink do?
#' @param X a [bigstatsr::FBM.code256] matrix (as found in the `genotypes`
#' slot of a [bigsnpr::bigSNP] object).
#' @param ind.row An optional vector of the row indices that are used.
#' If not specified, all rows are used. Don't use negative indices.
#' @param ind.col An optional vector of the column indices that are used. If not
#'  specified, all columns are used. Don't use negative indices.
#' @param as.counts whether the counts of similar alleles, rather than the proportion,
#' should be returned (FALSE by default). CURRENTLY ALWAYS RETURNS COUNTS!!!
#' @param block.size maximum number of columns read at once. Note that, to optimise the
#' speed of matrix operations, we have to store in memory 3 times the columns.
#' @returns a list of two [bigstatsr::FBM] matrices, one of counts of IBS by alleles (i.e. 2*n loci),
#' and one of valid alleles (i.e. 2*n_loci - 2*missing_loci)
#' @export

snp_ibs <- function(
  X,
  ind.row = bigstatsr::rows_along(X),
  ind.col = bigstatsr::cols_along(X),
  as.counts = FALSE,
  block.size = bigstatsr::block_size(nrow(X))
  ) {

  #check_args()

  n <- length(ind.row)
  # FBM matrix to count the IBS counts
  K <- bigstatsr::FBM(n, n, init = 0)
  # FBM matrix to store valid number of comparisons (i.e NOT NA)
  K2 <- bigstatsr::FBM(n, n, init = 0)
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
    ind.col.ind <- ind.col[ind] #subset given indices by the iteration indices
    increment_ibs_counts(K,
                         K2,
                         X_part_temp0,
                         X_part_temp1,
                         X_part_temp2,
                         X,
                         ind.row,
                         ind.col.ind)
  }

  return(list(ibs = K, valid_n = K2))
}


## convenience functions that are not exported by `bigstatsr`
CutBySize <- function (m, block.size, nb = ceiling(m/block.size))
{
  bigparallelr::split_len(m, nb_split = nb)
}

seq2<- function (lims)
{
  seq(lims[1], lims[2])
}


