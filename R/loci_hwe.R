#' Test Hardy-Weinberg equilibrium at each locus
#'
#' Return the p-value from an exact test of HWE.
#'
#' This function uses the original C++ algorithm from PLINK 1.90.
#'
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotypes` column of
#'   a [`gen_tibble`] object), or a [`gen_tibble`].
#' @param .col the column to be used when a tibble (or grouped tibble is passed
#' directly to the function). This defaults to "genotypes" and can only take
#' that value. There is no need for the user to set it, but it is included to
#' resolve certain tidyselect operations.
#' @param n_cores number of cores to be used, it defaults to
#'   [bigstatsr::nb_cores()]
#' @param block_size maximum number of loci read at once.
#' @param mid_p boolean on whether the mid-p value should be computed. Default
#'   is TRUE, as in PLINK.
#' @param type type of object to return, if using grouped method. One of "tidy",
#'   "list", or "matrix". Default is "tidy".
#'   Data can be grouped by more than one variable: groups are then labelled
#'   by the values of the grouping variables pasted with "." (e.g. `pop1.F`),
#'   and the "tidy" output has one column per grouping variable rather than a
#'   single `group` column. The `group` column, or the columns of the
#'   grouping variables, keep the type of the grouping variables (e.g. a
#'   factor keeps its levels), so that groups stay in their order when the
#'   output is sorted or plotted.
#' @param ... not used.
#' @returns a vector of probabilities from HWE exact test, one per locus. For
#'   a grouped [`gen_tibble`], empty groups (e.g. from
#'   `group_by(.drop = FALSE)`) are NA.
#' @author the C++ algorithm was written by Christopher Chang for PLINK 1.90,
#'   based on original code by Jan Wigginton (the code was released under GPL3).
#' @rdname loci_hwe
#' @export
#' @examples
#' \dontshow{
#' data.table::setDTthreads(2)
#' RhpcBLASctl::blas_set_num_threads(2)
#' RhpcBLASctl::omp_set_num_threads(2)
#' }
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # For HWE
#' example_gt %>% loci_hwe()
#'
#' # For loci_hwe per locus per population, use reframe
#' example_gt %>%
#'   group_by(population) %>%
#'   reframe(loci_hwe = loci_hwe(genotypes))
#'
loci_hwe <- function(.x, .col = "genotypes", ...) {
  UseMethod("loci_hwe", .x)
}


#' @export
#' @rdname loci_hwe
loci_hwe.tbl_df <- function(
  .x,
  .col = "genotypes",
  mid_p = TRUE,
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(nrow(.x), 1),
  ...
) {
  stopifnot_gen_tibble(.x)
  .col <- rlang::enquo(.col) %>%
    rlang::quo_get_expr() %>%
    rlang::as_string()
  # confirm that .col is "genotypes"
  if (.col != "genotypes") {
    stop("loci_hwe only works with the genotypes column")
  }
  loci_hwe(.x$genotypes,
    mid_p = mid_p, n_cores = n_cores,
    block_size = block_size, ...
  )
}


#' @export
#' @rdname loci_hwe
loci_hwe.vctrs_bigSNP <- function(
  .x,
  .col = "genotypes",
  mid_p = TRUE,
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(length(.x), 1),
  ...
) {
  rlang::check_dots_empty()
  stopifnot_diploid(.x)
  # rows (individuals) that we want to use
  rows_to_keep <- vctrs::vec_data(.x)
  # as long as we have more than one individual
  if (length(rows_to_keep) > 1) {
    # count genotypes and test HWE at each locus, treating all individuals as
    # a single group (the cpp code reads the FBM directly, in blocks of
    # block_size loci, using openMP for parallelisation)
    hwe_p <- gt_grouped_hwe(
      BM = attr(.x, "fbm"),
      rowInd = rows_to_keep,
      colInd = attr(.x, "loci")$big_index,
      groupIds = rep(0L, length(rows_to_keep)),
      ngroups = 1,
      midp = mid_p,
      ncores = n_cores,
      block_size = block_size
    )
    # drop the dimensions of the one column matrix (without copying it)
    dim(hwe_p) <- NULL
  } else if (length(rows_to_keep) == 0) {
    # no individuals (e.g. an empty group in group_map()), so no genotypes to
    # test, as for empty groups in the grouped method
    hwe_p <- rep(NA_real_, length(attr(.x, "loci")$big_index))
  } else {
    # if we have a single individual
    stop("Not implemented for a single individual")
  }
  hwe_p
}


#' @export
#' @rdname loci_hwe
loci_hwe.grouped_df <- function(
  .x,
  .col = "genotypes",
  mid_p = TRUE,
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(nrow(.x), 1), # nolint
  type = c("tidy", "list", "matrix"),
  ...
) {
  stopifnot_diploid(.x)

  .col <- rlang::enquo(.col) %>%
    rlang::quo_get_expr() %>%
    rlang::as_string()
  # confirm that .col is "genotypes"
  if (.col != "genotypes") {
    stop("loci_hwe only works with the genotypes column")
  }

  rlang::check_dots_empty()
  type <- match.arg(type)
  # count genotypes and test HWE at each locus for each group (the cpp code
  # reads the FBM directly, in blocks of block_size loci, and fills a single
  # matrix; openMP is used for parallelisation)
  hwe_mat <- gt_grouped_hwe( # nolint
    BM = .gt_get_fbm(.x),
    rowInd = .gt_fbm_rows(.x),
    colInd = .gt_fbm_cols(.x),
    groupIds = dplyr::group_indices(.x) - 1L,
    ngroups = dplyr::n_groups(.x),
    midp = mid_p,
    ncores = n_cores,
    block_size = block_size
  )

  groups <- .group_info(.x)
  # empty groups have no genotypes to test, so they are NA (rather than the
  # p-value of a test on zero counts)
  if (any(groups$sizes == 0)) {
    hwe_mat[, groups$sizes == 0] <- NA_real_
  }
  hwe_mat <- format_grouped_output(
    out_mat = hwe_mat,
    group_ids = groups$labels,
    loci_names = loci_names(.x),
    type = type,
    group_keys = groups$keys
  )
  return(hwe_mat)
}
