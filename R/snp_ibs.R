#' Compute the Identity by State Matrix for a bigSNP object
#'
#' This function computes the IBS matrix.
#'
#' Note that monomorphic sites are counted (as in PLINK). Remove monomorphic
#' sites beforehand if this is a concern.
#' @param X a [bigstatsr::FBM.code256] matrix (as found in the `genotypes`
#' slot of a [bigsnpr::bigSNP] object). Genotypes are decoded with its
#' `code256`, and values other than 0, 1 and 2 are treated as missing.
#' @param ind.row An optional vector of the row indices that are used. If not
#'   specified, all rows are used. Don't use negative indices.
#' @param ind.col An optional vector of the column indices that are used. If not
#'   specified, all columns are used. Don't use negative indices.
#' @param type one of:
#'   - "proportion" (the default): the proportion of alleles that are identical
#'   by state, a similarity (as `--distance ibs flat-missing` in PLINK, which
#'   writes it to a .mibs file).
#'   - "adjusted_counts": the allele count distance, a dissimilarity (as
#'   `--distance flat-missing` in PLINK, which writes it to a .dist file): the
#'   number of alleles that differ between two individuals, scaled up by the
#'   number of loci over the number of loci genotyped in both individuals, to
#'   correct for missing genotypes. It equals `2 * n_loci * (1 -
#'   proportion)`, where `n_loci` is the number of loci used, so it is 0 for
#'   identical individuals and at most `2 * n_loci`. Note that, without
#'   `flat-missing`, PLINK corrects for missing genotypes differently, which
#'   gives slightly different values when genotypes are missing.
#'   - "raw_counts": the counts of identical alleles and of non-missing
#'   alleles, from which the two other quantities are computed.
#' @param block.size maximum number of columns read at once. Genotypes are
#'   compared in a compact bit-packed form, so this has little effect on
#'   memory use.
#' @param n_cores number of cores to be used, it defaults to
#'   [bigstatsr::nb_cores()]
#' @returns for "proportion" and "adjusted_counts", a matrix with one row and
#'   one column per individual; for "raw_counts", a list of two
#'   [bigstatsr::FBM] matrices, `ibs` (the counts of alleles identical by
#'   state) and `valid_n` (the counts of non-missing alleles, i.e. twice the
#'   number of loci genotyped in both individuals).
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
  n_cores = bigstatsr::nb_cores()
) {
  # nolint end
  type <- match.arg(type)

  n <- length(ind.row)
  # no individuals (e.g. an empty group in group_map()): an empty matrix, as
  # FBMs can not be empty
  if (n == 0) {
    if (type == "raw_counts") {
      stop("type = 'raw_counts' needs at least one individual")
    }
    return(matrix(numeric(0), 0, 0))
  }
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
      # the allele count distance (as PLINK --distance flat-missing): the
      # number of alleles that differ, scaled up to all m loci to correct for
      # missing genotypes, i.e. 2 * m * (1 - ibs_prop)
      ibs_dist <- 2 * m * (1 - ibs_prop)
      return(ibs_dist)
    }
  }
}
