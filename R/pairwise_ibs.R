#' Compute the Identity by State Matrix for a `gen_tibble` object
#'
#' This function computes the IBS matrix.
#'
#' Note that monomorphic sites are currently considered. Remove monomorphic
#' sites before running pairwise_king if this is a concern.
#'
#' Imputed genotypes are used if `x` is set to use them (see
#' [gt_set_imputed()]), and treated as missing otherwise.
#' @param x a `gen_tibble` object.
#' @param as_matrix boolean, determining whether the results should be a square
#'   symmetrical matrix (TRUE), or a tidied tibble (FALSE, the default)
#' @param type one of "proportion" (equivalent to "ibs" in PLINK),
#'   "adjusted_counts" ("distance" in PLINK), and "raw_counts" (the counts of
#'   identical alleles and non-missing alleles, from which the two other
#'   quantities are computed)
#' @param block_size maximum number of loci read at once. Genotypes are
#'   compared in a compact bit-packed form, so this has little effect on
#'   memory use.
#' @param n_cores number of cores to be used (default 1).
#' @returns a [bigstatsr::FBM] of proportion or adjusted counts, or a list of
#'   two [bigstatsr::FBM] matrices, one of counts of IBS by alleles, and one of
#'   number of valid alleles (i.e. 2*n_loci - 2*missing_loci)
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' pairwise_ibs(example_gt, type = "proportion")
#'
#' # Alternatively, return a matrix
#' pairwise_ibs(example_gt, type = "proportion", as_matrix = TRUE)
#'
#' # Adjust block_size
#' pairwise_ibs(example_gt, block_size = 2)
#'
#' # Change type
#' pairwise_ibs(example_gt, type = "adjusted_counts")
#' pairwise_ibs(example_gt, type = "raw_counts")
pairwise_ibs <- function(
  x,
  as_matrix = FALSE,
  type = c(
    "proportion",
    "adjusted_counts",
    "raw_counts"
  ),
  block_size = bigstatsr::block_size(nrow(x)),
  n_cores = 1
) {
  stopifnot_diploid(x)
  type <- match.arg(type)
  X <- attr(x$genotypes, "fbm") # convenient pointer #nolint
  x_ind_col <- show_loci(x)$big_index
  x_ind_row <- vctrs::vec_data(x$genotypes)
  ibs_matrix <- snp_ibs(
    X,
    ind.row = x_ind_row,
    ind.col = x_ind_col,
    type = type,
    block.size = block_size,
    n_cores = n_cores
  )
  if (inherits(ibs_matrix, "matrix")) {
    dimnames(ibs_matrix) <- list(x$id, x$id)
  } else {
    # else if we have a list of two count matrices
    attr(ibs_matrix[[1]], "indiv_names") <- x$id
    attr(ibs_matrix[[2]], "indiv_names") <- x$id
  }

  if (type == "proportion" || type == "adjusted_counts") {
    if (as_matrix) {
      # add a class to indicate this is a pairwise matrix
      class(ibs_matrix) <- c("pairwise_matrix", class(ibs_matrix))
      return(ibs_matrix)
    } else {
      return(tidy_dist_matrix(ibs_matrix))
    }
  }

  ibs_matrix
}
