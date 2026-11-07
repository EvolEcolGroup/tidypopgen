#' Compute the Pairwise Allele Sharing Matrix for a bigSNP object
#'
#' This function computes the Allele Sharing matrix. Estimates Allele Sharing
#' (matching in `hierfstat`)) between pairs of individuals (for each locus,
#' gives 1 if the two individuals are homozygous for the same allele, 0 if they
#' are homozygous for a different allele, and 1/2 if at least one individual is
#' heterozygous. Matching is the average of these 0, 1/2 and 1s)
#'
#' @param X a [bigstatsr::FBM.code256] matrix (as found in the `genotypes`
#' slot of a [bigsnpr::bigSNP] object). Genotypes are decoded with its
#' `code256`, and values other than 0, 1 and 2 are treated as missing.
#' @param ind.row An optional vector of the row indices that are used. If not
#'   specified, all rows are used. Don't use negative indices.
#' @param ind.col An optional vector of the column indices that are used. If not
#'   specified, all columns are used. Don't use negative indices.
#' @param block.size maximum number of columns read at once. Genotypes are
#'   compared in a compact bit-packed form, so this has little effect on
#'   memory use.
#' @param n_cores number of cores to be used, it defaults to
#'   [bigstatsr::nb_cores()]
#' @returns a matrix of allele sharing between all pairs of individuals
#' @export
#' @seealso [pairwise_allele_sharing()] [hierfstat::matching()]
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' X <- attr(example_gt$genotypes, "fbm")
#' snp_allele_sharing(X)
#'
#' # Compute for individuals 1 to 5
#' snp_allele_sharing(X, ind.row = 1:5, ind.col = 1:5)
#'
#' # Adjust block size
#' snp_allele_sharing(X, block.size = 2)
#'
snp_allele_sharing <- function(
  X, # nolint start
  ind.row = bigstatsr::rows_along(X),
  ind.col = bigstatsr::cols_along(X),
  block.size = bigstatsr::block_size(nrow(X)),
  n_cores = bigstatsr::nb_cores()
) {
  # nolint end

  n <- length(ind.row)
  # no individuals (e.g. an empty group in group_map()): an empty matrix, as
  # FBMs can not be empty
  if (n == 0) {
    return(matrix(numeric(0), 0, 0))
  }
  # FBM matrix to store the tcrossproduct of (dos-1)
  dos_tcross <- bigstatsr::FBM(n, n, init = 0)
  # FBM matrix to store tcrossproduct of NA matrix
  # (0 is NA, 1 is a valid genotype)
  na_tcross <- bigstatsr::FBM(n, n, init = 0)
  m <- length(ind.col)

  intervals <- CutBySize(m, block.size)

  # for each block of loci, genotypes are packed into bits and compared
  # with popcounts in C++ (see src/pairwise_bitpacked.h), so no dense
  # copy of the block is needed
  for (j in bigstatsr::rows_along(intervals)) {
    ind <- seq2(intervals[j, ]) # this iteration indices
    ind.col.ind <- ind.col[ind] # subset given indices by the iteration indices #nolint
    increment_as_counts(
      dos_tcross,
      na_tcross,
      X,
      ind.row,
      ind.col.ind,
      ncores = n_cores
    )
  }
  # The allele sharing matrix is:
  # Mij <- 1/2 * (1+1/tcrossprod(na) * tcrossprod(dos-1)) #nolint
  # where dos is a dosage matrix (with missing values replaced by 1) and
  # na is a matrix with 0 for NAs and 1 for valid values
  # TODO this could be done in chunks to avoid bringing
  # everything to memory at once
  # (drop = FALSE keeps a 1 x 1 matrix for a single individual)
  denom <- na_tcross[, , drop = FALSE]
  num <- dos_tcross[, , drop = FALSE]
  res <- 0.5 * (1 + num / denom)
  res[denom == 0] <- NA_real_
  return(res)
}
