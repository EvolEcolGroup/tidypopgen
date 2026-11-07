#' Compute the KING-robust Matrix for a `gen_tibble` object
#'
#' This function computes the KING-robust estimator of kinship, reimplementing
#' the KING algorithm of Manichaikul et al. (2010).
#'
#' Imputed genotypes are used if `x` is set to use them (see
#' [gt_set_imputed()]), and treated as missing otherwise.
#'
#' `pairwise_king()` is an S3 generic. The method for a [`gen_tibble`]
#' computes the matrix with the method for its `genotypes` column (of class
#' `vctrs_bigSNP`), and labels it with the `id` of the individuals; grouping
#' is ignored. Other packages can add methods for their own classes of
#' genotypes (e.g. genotype likelihoods).
#'
#' @references Manichaikul, A. et al. (2010) Robust relationship inference in
#'   genome-wide association studies. Bioinformatics, 26(22), 2867–2873.
#'   https://doi.org/10.1093/bioinformatics/btq559.
#'
#' Note that monomorphic sites are currently considered. Remove monomorphic
#' sites before running pairwise_king if this is a concern.
#' @param x a `gen_tibble` object, or a vector of class `vctrs_bigSNP` (the
#'   `genotypes` column of a [`gen_tibble`], whose individuals are then
#'   labelled with the names of the vector).
#' @param as_matrix boolean, determining whether the results should be a square
#'   symmetrical matrix (TRUE), or a tidied tibble (FALSE, the default)
#' @param block_size maximum number of loci read at once. Genotypes are
#'   compared in a compact bit-packed form, so this has little effect on
#'   memory use.
#' @param n_cores number of cores to be used, it defaults to
#'   [bigstatsr::nb_cores()]
#' @param ... other arguments passed to specific methods, currently unused.
#' @returns a square symmetrical matrix of relationship coefficients between
#'   individuals if `as_matrix` is TRUE, or a tidied tibble of coefficients if
#'   `as_matrix` is FALSE.
#' @rdname pairwise_king
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # Compute the KING-robust matrix
#' pairwise_king(example_gt, as_matrix = TRUE)
#'
#' # Or return a tidy tibble
#' pairwise_king(example_gt, as_matrix = FALSE)
#'
#' # Adjust block_size
#' pairwise_king(example_gt, block_size = 2)
#'
pairwise_king <- function(x, ...) {
  UseMethod("pairwise_king", x)
}

#' @export
#' @rdname pairwise_king
pairwise_king.tbl_df <- function(
  x,
  as_matrix = FALSE,
  block_size = bigstatsr::block_size(nrow(x)),
  n_cores = bigstatsr::nb_cores(),
  ...
) {
  stopifnot_has_genotypes(x)
  # hand the genotypes column over to its method, which returns a
  # pairwise_matrix, and label it with the ids of the tibble (the matrix is
  # modified in place, without copies)
  king_matrix <- pairwise_king(
    x$genotypes,
    as_matrix = TRUE,
    block_size = block_size,
    n_cores = n_cores,
    ...
  )
  dimnames(king_matrix) <- list(x$id, x$id)
  if (as_matrix) {
    return(king_matrix)
  } else {
    class(king_matrix) <- NULL
    return(tidy_dist_matrix(king_matrix))
  }
}

#' @export
#' @rdname pairwise_king
pairwise_king.vctrs_bigSNP <- function(
  x,
  as_matrix = FALSE,
  block_size = bigstatsr::block_size(length(x)),
  n_cores = bigstatsr::nb_cores(),
  ...
) {
  rlang::check_dots_empty()
  stopifnot_diploid(x)
  king_matrix <- snp_king(
    attr(x, "fbm"),
    ind.row = vctrs::vec_data(x),
    ind.col = show_loci(x)$big_index,
    block.size = block_size,
    n_cores = n_cores
  )
  # the individuals are labelled with the names of the vector (the ids of the
  # gen_tibble)
  dimnames(king_matrix) <- list(names(x), names(x))
  if (as_matrix) {
    # add a class to indicate this is a pairwise matrix
    class(king_matrix) <- c("pairwise_matrix", class(king_matrix))
    return(king_matrix)
  } else {
    return(tidy_dist_matrix(king_matrix))
  }
}


# KING-robust kinship within each group of a grouped gen_tibble, as a list of
# matrices (one per group, in the order of the groups, with NULL for empty
# groups), identical to the output of pairwise_king(as_matrix = TRUE) on each
# group. Rather than reading all the loci once for each group, the genotypes
# of all individuals are packed once, and only the pairs of individuals
# within each group are compared (see king_counts_within_groups()).
.pairwise_king_by_group <- function(
  .x,
  block_size = bigstatsr::block_size(nrow(.x)),
  n_cores = bigstatsr::nb_cores()
) {
  stopifnot_diploid(.x)
  X <- attr(.x$genotypes, "fbm") # nolint
  # the rows of each group, and the individuals (rows of the FBM) ordered by
  # group (keeping their order within each group)
  group_rows <- dplyr::group_rows(.x)
  sizes <- lengths(group_rows)
  king <- vector("list", length(sizes))
  if (sum(sizes) == 0) {
    return(king)
  }
  ind_row <- vctrs::vec_data(.x$genotypes)[unlist(group_rows)]
  # the position (in ind_row) of the last individual of the group of each
  # individual
  group_end <- rep(cumsum(sizes), sizes)
  ind_col <- show_loci(.x)$big_index
  # counts for each pair within each group, summed over blocks of loci (the
  # counts are integers, so the sums are exact)
  counts <- 0
  intervals <- CutBySize(length(ind_col), block_size)
  for (j in bigstatsr::rows_along(intervals)) {
    counts <- counts + king_counts_within_groups(X,
      rowInd = ind_row,
      colInd = ind_col[seq2(intervals[j, ])],
      groupEnd = group_end,
      ncores = n_cores
    )
  }
  # the pairs of each group (i <= j) are consecutive columns of counts
  n_pairs <- sizes * (sizes + 1) / 2
  pair_offset <- c(0, cumsum(n_pairs))
  for (g in which(sizes > 0)) {
    k <- sizes[g]
    pairs <- pair_offset[g] + seq_len(n_pairs[g])
    i <- rep(seq_len(k), times = k:1)
    j <- i + sequence(k:1) - 1L
    # the numerator of KING (symmetric), and the number of loci at which the
    # individual of the row is heterozygous and the one of the column valid
    numerator <- matrix(0, k, k)
    numerator[cbind(i, j)] <- numerator[cbind(j, i)] <-
      counts[1, pairs] - 2 * counts[2, pairs]
    n_het <- matrix(0, k, k)
    n_het[cbind(i, j)] <- counts[3, pairs]
    n_het[cbind(j, i)] <- counts[4, pairs]
    n_het_t <- t(n_het)
    # as in snp_king()
    king_g <- numerator /
      (2 * pmin(n_het, n_het_t)) +
      0.5 -
      0.25 * (n_het + n_het_t) / pmin(n_het, n_het_t)
    ids <- .x$id[group_rows[[g]]]
    dimnames(king_g) <- list(ids, ids)
    class(king_g) <- c("pairwise_matrix", class(king_g))
    king[[g]] <- king_g
  }
  king
}
