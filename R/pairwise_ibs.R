#' Compute the Identity by State Matrix for a `gen_tibble` object
#'
#' This function computes the IBS matrix.
#'
#' Note that monomorphic sites are currently considered (as in PLINK). Remove
#' monomorphic sites before running `pairwise_ibs()` if this is a concern.
#'
#' Imputed genotypes are used if `x` is set to use them (see
#' [gt_set_imputed()]), and treated as missing otherwise.
#'
#' Missing genotypes are handled as in PLINK with the `flat-missing` modifier
#' (`--distance ibs flat-missing` for "proportion", `--distance
#' flat-missing` for "adjusted_counts"): each pair of individuals is compared
#' only on the loci genotyped in both, and all loci have the same weight. By
#' default (without `flat-missing`), PLINK 1.9 instead divides each distance
#' by one minus the sum of the average contributions to the distance of the
#' loci missing in the pair, which depend on their allele frequencies. This
#' default correction is not available in `pairwise_ibs()`, so, when
#' genotypes are missing, its values are only comparable to those of PLINK
#' run with `flat-missing`. Without missing genotypes, both corrections give
#' the same values.
#' @param x a `gen_tibble` object.
#' @param as_matrix boolean, determining whether the results should be a square
#'   symmetrical matrix (TRUE), or a tidied tibble (FALSE, the default)
#' @param type one of:
#'   - "proportion" (the default): the proportion of alleles that are identical
#'   by state, a similarity (as `--distance ibs flat-missing` in PLINK, which
#'   writes it to a .mibs file).
#'   - "adjusted_counts": the allele count distance, a dissimilarity (as
#'   `--distance flat-missing` in PLINK, which writes it to a .dist file): the
#'   number of alleles that differ between two individuals, scaled up by the
#'   number of loci over the number of loci genotyped in both individuals, to
#'   correct for missing genotypes. It equals `2 * n_loci * (1 -
#'   proportion)`, so it is 0 for identical individuals and at most `2 *
#'   n_loci`. PLINK's default correction for missing genotypes (without
#'   `flat-missing`) is not available (see Details).
#'   - "raw_counts": the counts of identical alleles and of non-missing
#'   alleles, from which the two other quantities are computed.
#' @param block_size maximum number of loci read at once. Genotypes are
#'   compared in a compact bit-packed form, so this has little effect on
#'   memory use.
#' @param n_cores number of cores to be used, it defaults to
#'   [bigstatsr::nb_cores()]
#' @returns for "proportion" and "adjusted_counts", a tidy tibble with one row
#'   per pair of individuals (columns `item1`, `item2` and `value`), or a
#'   square matrix if `as_matrix = TRUE`; for "raw_counts", a list of two
#'   [bigstatsr::FBM] matrices, `ibs` (the counts of alleles identical by
#'   state) and `valid_n` (the counts of non-missing alleles, i.e. twice the
#'   number of loci genotyped in both individuals).
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
  n_cores = bigstatsr::nb_cores()
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
