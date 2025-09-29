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
#' "list", or "matrix". Default is "tidy".
#' @param ... not used.
#' @returns a vector of probabilities from HWE exact test, one per locus
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
loci_hwe.tbl_df <- function(.x, .col = "genotypes", mid_p = TRUE, ...) {
  stopifnot_gen_tibble(.x)
  .col <- rlang::enquo(.col) %>%
    rlang::quo_get_expr() %>%
    rlang::as_string()
  # confirm that .col is "genotypes"
  if (.col != "genotypes") {
    stop("loci_hwe only works with the genotypes column")
  }
  loci_hwe(.x$genotypes, mid_p = mid_p, ...)
}


#' @export
#' @rdname loci_hwe
loci_hwe.vctrs_bigSNP <- function(.x, .col = "genotypes", mid_p = TRUE, ...) {
  rlang::check_dots_empty()
  stopifnot_diploid(.x)
  # get the FBM
  geno_fbm <- attr(.x, "fbm")
  # rows (individuals) that we want to use
  rows_to_keep <- vctrs::vec_data(.x)
  # as long as we have more than one individual
  if (length(rows_to_keep) > 1) {
    # col hwe for submatrix (some rows, and some columns)
    col_hwe_sub <- function(X, ind, rows_to_keep) { # nolint
      # apply(X[rows_to_keep, ind], 2, HWExact_geno_vec) #nolint
      geno_counts <- bigstatsr::big_counts(
        X,
        ind.row = rows_to_keep,
        ind.col = ind
      )
      hwe_on_matrix(geno_counts = geno_counts, midp = mid_p)
    }
    hwe_p <- bigstatsr::big_apply(
      geno_fbm,
      a.FUN = col_hwe_sub,
      rows_to_keep = rows_to_keep,
      ind = attr(.x, "loci")$big_index,
      a.combine = "c"
    )
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
    ...) {
  stopifnot_diploid(.x)

  .col <- rlang::enquo(.col) %>%
    rlang::quo_get_expr() %>%
    rlang::as_string()
  # confirm that .col is "genotypes"
  if (.col != "genotypes") {
    stop("loci_hwe only works with the genotypes column")
  }

  # check that we only have one grouping variable
  if (length(.x %>% dplyr::group_vars()) > 1) {
    stop("loci_hwe only works with one grouping variable")
  }
  rlang::check_dots_empty()
  type <- match.arg(type)
  geno_fbm <- .gt_get_fbm(.x)
  rows_to_keep <- .gt_fbm_rows(.x)
  hwe_p_sub <- function(geno_fbm, ind, rows_to_keep) {
    gt_grouped_hwe( # nolint
      BM = geno_fbm,
      rowInd = rows_to_keep,
      colInd = ind,
      groupIds = dplyr::group_indices(.x) - 1,
      ngroups = max(dplyr::group_indices(.x)),
      midp = mid_p
    )
  }

  hwe_mat <- bigstatsr::big_apply(
    geno_fbm,
    a.FUN = hwe_p_sub,
    rows_to_keep = rows_to_keep,
    ind = attr(.x$genotypes, "loci")$big_index,
    ncores = 1, # parallelisation is used within the function
    block.size = block_size,
    a.combine = "rbind"
  )

  hwe_mat <- format_grouped_output(
    out_mat = hwe_mat,
    group_ids = dplyr::group_keys(.x) %>% pull(1),
    loci_names = loci_names(.x),
    type = type
  )
  return(hwe_mat)
}
