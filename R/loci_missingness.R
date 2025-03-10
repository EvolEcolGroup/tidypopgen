#' Estimate missingness at each locus
#'
#' Estimate the rate of missingness at each locus.
#'
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotypes` column of
#'   a [`gen_tibble`] object), or a [`gen_tibble`].
#' @param as_counts boolean defining whether the count of NAs (rather than the
#'   rate) should be returned. It defaults to FALSE (i.e. rates are returned by
#'   default).
#' @param n_cores number of cores to be used, it defaults to
#'   [bigstatsr::nb_cores()]
#' @param block_size maximum number of loci read at once.
#' @param ... other arguments passed to specific methods.
#' @returns a vector of frequencies, one per locus
#' @rdname loci_missingness
#' @export
loci_missingness <- function(.x, as_counts = FALSE, block_size, ...) {
  UseMethod("loci_missingness", .x)
}

#' @export
#' @rdname loci_missingness
loci_missingness.tbl_df <- function(
    .x,
    as_counts = FALSE,
    # the bigapply that splits in blocks is not
    # multithreaded, as we use the multiple
    # threads for openMP,
    block_size = bigstatsr::block_size(nrow(attr(.x$genotypes, "loci")), 1), # nolint
    ...) {
  # TODO this is a hack to deal with the class being dropped when going through
  # group_map
  stopifnot_gen_tibble(.x)
  loci_missingness(
    .x$genotypes,
    as_counts = as_counts,
    block_size = block_size,
    ...
  )
}


#' @export
#' @rdname loci_missingness
loci_missingness.vctrs_bigSNP <- function(
    .x,
    as_counts = FALSE,
    block_size = bigstatsr::block_size(nrow(attr(.x, "loci")), 1), # nolint
    ...) {
  rlang::check_dots_empty()
  # get the FBM
  geno_fbm <- attr(.x, "bigsnp")$genotypes
  # rows (individuals) that we want to use
  rows_to_keep <- vctrs::vec_data(.x)
  # as long as we have more than one individual
  if (length(rows_to_keep) > 1) {
    # internal function that can be used with a big_apply
    count_na_sub <- function(BM, ind, rows_to_keep) {
      # nolint
      n_na <- bigstatsr::big_counts(BM, ind.row = rows_to_keep, ind.col = ind)
      n_na <- n_na[nrow(n_na), ] # this should work also with polyploids
      n_na
    }
    n_na <- bigstatsr::big_apply(
      geno_fbm,
      a.FUN = count_na_sub,
      rows_to_keep = rows_to_keep,
      ind = attr(.x, "loci")$big_index,
      ncores = 1, # parallelisation is used within the function
      block.size = block_size,
      a.combine = "c"
    )
    if (!as_counts) {
      n_na <- n_na / length(rows_to_keep)
    }
  } else {
    # if we have a single individual
    n_na <- geno_fbm[rows_to_keep, attr(.x, "loci")$big_index]
  }
  n_na
}

#' @export
#' @rdname loci_missingness
loci_missingness.grouped_df <- function(
    .x,
    as_counts = FALSE,
    block_size = bigstatsr::block_size(nrow(attr(.x, "loci")), 1), # nolint
    n_cores = bigstatsr::nb_cores(),
    ...) {
  rlang::check_dots_empty()
  geno_fbm <- .gt_get_bigsnp(.x)$genotypes
  rows_to_keep <- .gt_bigsnp_rows(.x)
  count_na_sub <- function(geno_fbm, ind, rows_to_keep) {
    na_mat <- gt_grouped_missingness(
      # nolint
      BM = geno_fbm,
      rowInd = rows_to_keep,
      colInd = ind,
      groupIds = dplyr::group_indices(.x) - 1,
      ngroups = max(dplyr::group_indices(.x)),
      ncores = n_cores
    )
  }

  na_mat <- bigstatsr::big_apply(
    geno_fbm,
    a.FUN = count_na_sub,
    rows_to_keep = rows_to_keep,
    ind = attr(.x$genotypes, "loci")$big_index,
    ncores = 1, # parallelisation is used within the function
    block.size = block_size,
    a.combine = "rbind"
  )

  group_sizes <- tally(.x) %>% dplyr::pull(dplyr::all_of("n"))
  if (!as_counts) {
    na_mat <- sweep(na_mat, MARGIN = 2, STATS = group_sizes, FUN = "/")
  }

  # return a list to mimic a group_map
  lapply(seq_len(ncol(na_mat)), function(i) na_mat[, i])
}
