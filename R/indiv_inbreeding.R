#' Individual inbreeding coefficient
#'
#' This function calculates the inbreeding coefficient for each individual based
#' on the beta estimate from Weir and Goudet (2017).
#'
#' `indiv_inbreeding()` is an S3 generic, with methods for [`gen_tibble`]
#' objects (grouped or not) and for their `genotypes` column. Other packages
#' can add methods for their own classes (e.g. tibbles of genotype
#' likelihoods).
#'
#' For a grouped [`gen_tibble`], the coefficients are computed within each
#' group, giving the same results as
#' `group_map(.x, ~ indiv_inbreeding(.x))`, but much faster: the genotypes are
#' read only once, and only the pairs of individuals within each group are
#' compared.
#'
#' @references Weir, BS and Goudet J (2017) A Unified Characterization of
#'   Population Structure and Relatedness. Genetics (2017) 206:2085
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotype` column of
#'   a [`gen_tibble`] object), or a [`gen_tibble`]. The method for the
#'   `genotypes` column needs `allele_sharing_mat`.
#' @param method currently only "WG17" (for Weir and Goudet 2017).
#' @param allele_sharing_mat optional and only relevant for "WG17", the matrix
#'   of Allele Sharing returned by [pairwise_allele_sharing()] with
#'   `as_matrix=TRUE`. As a number of statistics can be derived from the Allele
#'   Sharing matrix, it is sometimes more efficient to pre-compute this matrix.
#'   It is not possible to use this with grouped tibbles.
#' @param ... other arguments passed to specific methods, currently unused.
#' @return a numeric vector of inbreeding coefficients. For a grouped
#'   [`gen_tibble`], a list of such vectors, one per group in the order of
#'   [dplyr::group_keys()] (the coefficients are computed within each group;
#'   an empty group gives an empty vector).
#' @rdname indiv_inbreeding
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' example_gt %>% indiv_inbreeding(method = "WG17")
#'
indiv_inbreeding <- function(.x, ...) {
  UseMethod("indiv_inbreeding", .x)
}

#' @export
#' @rdname indiv_inbreeding
indiv_inbreeding.gen_tbl <- function(
  .x,
  method = c("WG17"),
  allele_sharing_mat = NULL,
  ...
) {
  stopifnot_has_genotypes(.x)
  # the matrix of allele sharing is computed from the whole gen_tibble (it is
  # labelled with its ids), then the genotypes column is handed over to its
  # method
  if (is.null(allele_sharing_mat)) {
    allele_sharing_mat <- pairwise_allele_sharing(.x, as_matrix = TRUE)
  }
  indiv_inbreeding(
    .x$genotypes,
    method = method,
    allele_sharing_mat = allele_sharing_mat,
    ...
  )
}

#' @export
#' @rdname indiv_inbreeding
indiv_inbreeding.vctrs_bigSNP <- function(
  .x,
  method = c("WG17"),
  allele_sharing_mat = NULL,
  ...
) {
  rlang::check_dots_empty()
  stopifnot_diploid(.x)
  if (is.null(allele_sharing_mat)) {
    # TODO deal with the fact that allele sharing can only be computed
    # for the gen_tibble
    # TODO modify allele sharing code to be able to create a matrix just
    # from the genotype column
    stop(paste(
      "this method has to be applied to a gen_tibble, not just",
      "the genotype column"
    ))
  }

  # taken from beta.dosage in hierfstat
  mii <- diag(allele_sharing_mat) # diagonal of the matrix
  diag(allele_sharing_mat) <- NA
  mb <- mean(allele_sharing_mat, na.rm = TRUE)
  indiv_inb <- ((mii * 2 - 1) - mb) / (1 - mb)
  return(indiv_inb)
}

#' @export
#' @rdname indiv_inbreeding
indiv_inbreeding.grouped_gen_tbl <- function(
  .x,
  method = c("WG17"),
  allele_sharing_mat = NULL,
  ...
) {
  rlang::check_dots_empty()
  stopifnot_diploid(.x)
  if (!is.null(allele_sharing_mat)) {
    stop("allele_sharing_mat can not be provided for grouped_df objects")
  }
  method <- match.arg(method)
  .indiv_inbreeding_by_group(.x)
}

# WG17 inbreeding within each group of a grouped gen_tibble, as a list of
# vectors (one per group, in the order of the groups), identical to the
# output of indiv_inbreeding() on each group (an empty vector for empty
# groups). Rather than computing a matrix of allele sharing for each group,
# the genotypes of all individuals are packed once, and only the pairs of
# individuals within each group are compared (see as_counts_within_groups()).
.indiv_inbreeding_by_group <- function(
  .x,
  block_size = bigstatsr::block_size(nrow(.x)),
  n_cores = bigstatsr::nb_cores()
) {
  X <- attr(.x$genotypes, "fbm") # nolint
  # the rows of each group, and the individuals (rows of the FBM) ordered by
  # group (keeping their order within each group)
  group_rows <- dplyr::group_rows(.x)
  sizes <- lengths(group_rows)
  # empty groups have an empty vector, as indiv_inbreeding() on a gen_tibble
  # without individuals
  inbreeding <- rep(list(numeric(0)), length(sizes))
  if (sum(sizes) == 0) {
    return(inbreeding)
  }
  ind_row <- vctrs::vec_data(.x$genotypes)[unlist(group_rows)]
  # the position (in ind_row) of the last individual of the group of each
  # individual
  group_end <- rep(cumsum(sizes), sizes)
  ind_col <- show_loci(.x$genotypes)$big_index
  # counts for each pair within each group, summed over blocks of loci (the
  # counts are integers, so the sums are exact)
  counts <- 0
  intervals <- CutBySize(length(ind_col), block_size)
  for (j in bigstatsr::rows_along(intervals)) {
    counts <- counts + as_counts_within_groups(X,
      rowInd = ind_row,
      colInd = ind_col[seq2(intervals[j, ])],
      groupEnd = group_end,
      ncores = n_cores
    )
  }
  # allele sharing of each pair, as in snp_allele_sharing() (NA for pairs
  # without loci valid in both)
  ashare <- 0.5 * (1 + counts[1, ] / counts[2, ])
  ashare[counts[2, ] == 0] <- NA_real_
  # the pairs of each group (i <= j) are consecutive columns of counts
  n_pairs <- sizes * (sizes + 1) / 2
  pair_offset <- c(0, cumsum(n_pairs))
  for (g in which(sizes > 0)) {
    k <- sizes[g]
    ashare_g <- ashare[pair_offset[g] + seq_len(n_pairs[g])]
    # whether each pair (in the order of the packed triangle) is on the
    # diagonal, i.e. an individual with itself
    on_diag <- sequence(k:1) == 1L
    # as in the method for the genotypes column (from beta.dosage in
    # hierfstat): the matrix is symmetric, so the mean of its off-diagonal
    # elements is the mean of one of its triangles (NaN for a single
    # individual)
    mii <- ashare_g[on_diag]
    mb <- mean(ashare_g[!on_diag], na.rm = TRUE)
    inbreeding_g <- ((mii * 2 - 1) - mb) / (1 - mb)
    names(inbreeding_g) <- .x$id[group_rows[[g]]]
    inbreeding[[g]] <- inbreeding_g
  }
  inbreeding
}

# a grouped data frame that has lost the class of its tibble (e.g. a grouped
# gen_tibble returned by summarise()): restore the class and dispatch again
#' @export
#' @noRd
indiv_inbreeding.grouped_df <- function(.x, ...) {
  indiv_inbreeding(restore_genotypes_tbl(.x, "indiv_inbreeding"), ...)
}

# a data frame that has lost the class of its tibble (e.g. a gen_tibble
# returned by pick() or reframe()): restore the class and dispatch again
#' @export
#' @noRd
indiv_inbreeding.data.frame <- function(.x, ...) {
  indiv_inbreeding(restore_genotypes_tbl(.x, "indiv_inbreeding"), ...)
}
