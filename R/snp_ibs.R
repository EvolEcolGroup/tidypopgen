#' Compute the Identity by State Matrix for a bigSNP object
#'
#' This function computes the IBS matrix.
#'
#' Note that monomorphic sites are currently counted. Should we filter them
#' beforehand? What does plink do?
#' @param X a [bigstatsr::FBM.code256] matrix (as found in the `genotypes`
#' slot of a [bigsnpr::bigSNP] object).
#' @param ind.row An optional vector of the row indices that are used. If not
#'   specified, all rows are used. Don't use negative indices.
#' @param ind.col An optional vector of the column indices that are used. If not
#'   specified, all columns are used. Don't use negative indices.
#' @param type one of "proportion" (equivalent to "ibs" in PLINK),
#'   "adjusted_counts" ("distance" in PLINK), and "raw_counts" (the counts of
#'   identical alleles and non-missing alleles, from which the two other
#'   quantities are computed)
#' @param block.size maximum number of columns read at once. Genotypes are
#'   compared in a compact bit-packed form, so this has little effect on
#'   memory use.
#' @param n_cores number of cores to be used (default 1).
#' @returns if as.counts = TRUE function returns a list of two [bigstatsr::FBM]
#'   matrices, one of counts of IBS by alleles (i.e. 2*n loci), and one of valid
#'   alleles (i.e. 2 * n_loci - 2 * missing_loci). If as.counts = FALSE returns
#'   a single matrix of IBS proportions.
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' X <- attr(example_gt$genotypes, "fbm")
#' snp_ibs(X)
#'
#' # Compute for individuals 1 to 5
#' snp_ibs(X, ind.row = 1:5, ind.col = 1:5)
#'
#' # Adjust block.size
#' snp_ibs(X, block.size = 2)
#'
#' # Change type
#' snp_ibs(X, type = "proportion")
#' snp_ibs(X, type = "adjusted_counts")
#' snp_ibs(X, type = "raw_counts")
#'
snp_ibs <- function(
  X, # nolint start
  ind.row = bigstatsr::rows_along(X),
  ind.col = bigstatsr::cols_along(X),
  type = c("proportion", "adjusted_counts", "raw_counts"),
  block.size = bigstatsr::block_size(nrow(X)),
  n_cores = 1
) {
  # nolint end
  type <- match.arg(type)

  n <- length(ind.row)
  # FBM matrix to count the IBS counts
  IBS <- bigstatsr::FBM(n, n, init = 0) # nolint
  # FBM matrix to store valid number of comparisons (i.e NOT NA)
  IBS_valid_loci <- bigstatsr::FBM(n, n, init = 0) # nolint
  m <- length(ind.col)

  intervals <- CutBySize(m, block.size)

  # for each block of loci, genotypes are packed into bits and compared
  # with popcounts in C++ (see src/pairwise_bitpacked.h), so no dense
  # copy of the block is needed
  for (j in bigstatsr::rows_along(intervals)) {
    ind <- seq2(intervals[j, ]) # this iteration indices
    ind.col.ind <- ind.col[ind] # subset given indices by the iteration indices #nolint
    increment_ibs_counts(
      IBS,
      IBS_valid_loci,
      X,
      ind.row,
      ind.col.ind,
      ncores = n_cores
    )
  }

  if (type == "raw_counts") {
    return(list(ibs = IBS, valid_n = IBS_valid_loci))
  } else {
    # IBS proportion
    # get the means of each column
    divide_sub <- function(X, ind, Y) (X[, ind] / Y[, ind]) # nolint
    ibs_prop <- bigstatsr::big_apply(
      IBS,
      a.FUN = divide_sub,
      Y = IBS_valid_loci,
      a.combine = "cbind"
    )
    if (type == "proportion") {
      return(ibs_prop)
    } else {
      # for adjusted counts
      ibs_adj <- ibs_prop * m
      return(ibs_adj)
    }
  }
}
