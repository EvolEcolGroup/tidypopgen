#' Estimate nucleotide diversity (pi) at each locus
#'
#' Estimate nucleotide diversity (pi) at each locus, accounting for missing
#' values. This uses the formula: c_0 * c_1 / (n * (n-1) / 2)
#'
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotypes` column of
#'   a [`gen_tibble`] object), or a [`gen_tibble`].
#' @param n_cores number of cores to be used, it defaults to
#'   [bigstatsr::nb_cores()]
#' @param block_size maximum number of loci read at once.
#' @param ... other arguments passed to specific methods, currently unused.
#' @returns a vector of frequencies, one per locus
#' @rdname loci_pi
#' @export
loci_pi <- function(.x, n_cores, block_size, ...) {
  UseMethod("loci_pi", .x)
}

#' @export
#' @rdname loci_pi
loci_pi.tbl_df <- function(
    .x,
    # multicore is used by openMP within the
    # freq cpp function
    n_cores = bigstatsr::nb_cores(),
    # the bigapply that splits in blocks is not
    # multithreaded, as we use the multiple threads
    # for openMP
    block_size = bigstatsr::block_size(nrow(.x), 1),
    ...) {
  # TODO this is a hack to deal with the class being dropped when going
  # through group_map
  stopifnot_gen_tibble(.x)
  loci_pi(.x$genotypes)
}


#' @export
#' @rdname loci_pi
loci_pi.vctrs_bigSNP <- function(
    .x,
    n_cores = bigstatsr::nb_cores(),
    block_size = bigstatsr::block_size(length(.x), 1),
    ...) {
  rlang::check_dots_empty()

  stopifnot_diploid(.x)
  # if we have diploid
  # get the FBM
  geno_fbm <- attr(.x, "bigsnp")$genotypes
  # rows (individuals) that we want to use
  rows_to_keep <- vctrs::vec_data(.x)
  # as long as we have more than one individual
  if (length(rows_to_keep) > 1) {
    # internal function that can be used with a big_apply #nolint start
    gt_pi_sub <- function(BM, ind, rows_to_keep) {
      gt_pi_diploid(
        BM = BM,
        rowInd = rows_to_keep,
        colInd = ind,
        ncores = n_cores
      )
    } # nolint end
    pi <- bigstatsr::big_apply(
      geno_fbm,
      a.FUN = gt_pi_sub,
      rows_to_keep = rows_to_keep,
      ind = attr(.x, "loci")$big_index,
      ncores = 1, # parallelisation is used within the function
      block.size = block_size,
      a.combine = "c"
    )
  } else {
    # if we have a single individual
    # pi does not really make sense for a single individual
    pi <- rep(NA_real_, length(attr(.x, "loci")$big_index))
  }
  pi
}

#' @export
#' @rdname loci_pi
loci_pi.grouped_df <- function(
    .x,
    n_cores = bigstatsr::nb_cores(),
    block_size = bigstatsr::block_size(nrow(.x), 1),
    ...) {
  rlang::check_dots_empty()
  stopifnot_diploid(.x)
  geno_fbm <- .gt_get_bigsnp(.x)$genotypes
  # rows (individuals) that we want to use
  rows_to_keep <- vctrs::vec_data(.x$genotypes)

  # internal function that can be used with a big_apply #nolint start
  gt_group_pi_sub <- function(BM, ind, rows_to_keep) {
    freq_mat <- gt_grouped_pi_diploid(
      BM = BM,
      rowInd = rows_to_keep,
      colInd = ind,
      groupIds = dplyr::group_indices(.x) - 1,
      ngroups = max(dplyr::group_indices(.x)),
      ncores = n_cores
    )$pi
  } # nolint end
  pi_mat <- bigstatsr::big_apply(
    geno_fbm,
    a.FUN = gt_group_pi_sub,
    rows_to_keep = rows_to_keep,
    ind = show_loci(.x)$big_index,
    ncores = 1, # we only use 1 cpu, we let openMP use multiple cores
    # in the cpp code
    block.size = block_size,
    a.combine = "rbind"
  )
  # return a list to mimic a group_map
  lapply(seq_len(ncol(pi_mat)), function(i) pi_mat[, i])
}
