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
#' group (the beta estimate is relative to the other individuals of the
#' group), giving the same values as
#' `group_map(.x, ~ indiv_inbreeding(.x))`, but much faster: the genotypes are
#' read only once, and only the pairs of individuals within each group are
#' compared. As for an ungrouped [`gen_tibble`], the result has one value per
#' individual, in the order of the rows, so that it can be added as a column
#' (e.g. `grouped_gt$inbreeding <- indiv_inbreeding(grouped_gt)`).
#'
#' The method for the `genotypes` column computes the matrix of allele
#' sharing of the individuals in the column, so that `indiv_inbreeding()` can
#' be used within `dplyr` verbs, e.g.
#' `mutate(inbreeding = indiv_inbreeding(genotypes))`. Within a grouped
#' `mutate()`, `dplyr` calls it on the individuals of each group in turn, so
#' it gives the same values as the method for a grouped [`gen_tibble`], but
#' it reads the genotypes once for each group.
#'
#' @references Weir, BS and Goudet J (2017) A Unified Characterization of
#'   Population Structure and Relatedness. Genetics (2017) 206:2085
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotype` column of
#'   a [`gen_tibble`] object), or a [`gen_tibble`].
#' @param method currently only "WG17" (for Weir and Goudet 2017).
#' @param allele_sharing_mat optional and only relevant for "WG17", the matrix
#'   of Allele Sharing returned by [pairwise_allele_sharing()] with
#'   `as_matrix=TRUE`. As a number of statistics can be derived from the Allele
#'   Sharing matrix, it is sometimes more efficient to pre-compute this matrix.
#'   It is not possible to use this with grouped tibbles.
#' @param type only for a grouped [`gen_tibble`]: "vector" (the default), one
#'   value per individual, in the order of the rows; or "list" (deprecated),
#'   one vector per group, in the order of [dplyr::group_keys()].
#' @param ... other arguments passed to specific methods, currently unused.
#' @return a numeric vector of inbreeding coefficients, one per individual in
#'   the order of the rows, named by the ids of the individuals (when
#'   `allele_sharing_mat` is given, by its row names). For a grouped
#'   [`gen_tibble`], the coefficients are computed within each group (with
#'   `type = "list"`, deprecated, they are returned as a list of vectors, one
#'   per group; an empty group gives an empty vector).
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
    # the matrix of allele sharing of the individuals in this column (e.g. of
    # one group, within a grouped mutate()), labelled with their ids
    allele_sharing_mat <- pairwise_allele_sharing(.x, as_matrix = TRUE)
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
  type = c("vector", "list"),
  ...
) {
  rlang::check_dots_empty()
  stopifnot_diploid(.x)
  if (!is.null(allele_sharing_mat)) {
    stop("allele_sharing_mat can not be provided for grouped_df objects")
  }
  method <- match.arg(method)
  type <- match.arg(type)
  # the coefficients within each group, one vector per group (in the order of
  # the groups), each named by the ids of the individuals of its group
  inbreeding <- .indiv_inbreeding_by_group(.x)
  if (type == "list") {
    # the output of tidypopgen 1.0.0, kept for compatibility
    lifecycle::deprecate_soft(
      "1.1.0", I('`indiv_inbreeding(type = "list")`'),
      details = paste(
        "Use the default `type = \"vector\"`, which gives one value per",
        "individual in the order of the rows."
      )
    )
    return(inbreeding)
  }
  # back to the order of the rows: the values of each group go to the rows
  # of its individuals (every row belongs to exactly one group)
  inbreeding_rows <- numeric(nrow(.x))
  inbreeding_rows[unlist(dplyr::group_rows(.x))] <-
    unlist(inbreeding, use.names = FALSE)
  names(inbreeding_rows) <- .x$id
  inbreeding_rows
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
