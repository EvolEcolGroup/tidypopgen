#' Estimate missingness at each locus
#'
#' Estimate the rate of missingness at each locus. This function has an
#' efficient method to support grouped `gen_tibble` objects, which can return a
#' tidied tibble, a list, or a matrix.
#'
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotypes` column of
#'   a [`gen_tibble`] object), or a [`gen_tibble`].
#' @param .col the column to be used when a tibble (or grouped tibble is passed
#' directly to the function). This defaults to "genotypes" and can only take
#' that value. There is no need for the user to set it, but it is included to
#' resolve certain tidyselect operations.
#' @param as_counts boolean defining whether the count of NAs (rather than the
#'   rate) should be returned. It defaults to FALSE (i.e. rates are returned by
#'   default).
#' @param n_cores number of cores to be used, it defaults to
#'   [bigstatsr::nb_cores()]
#' @param block_size maximum number of loci read at once.
#' @param type type of object to return, if using grouped method. One of "tidy",
#' "list", or "matrix". Default is "tidy".
#' @param ... other arguments passed to specific methods.
#' @returns a vector of frequencies, one per locus
#' @rdname loci_missingness
#' @export
#' @examplesIf all(rlang::is_installed(c("RhpcBLASctl", "data.table")))
#' \dontshow{
#' data.table::setDTthreads(2)
#' RhpcBLASctl::blas_set_num_threads(2)
#' RhpcBLASctl::omp_set_num_threads(2)
#' }
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # For missingness
#' example_gt %>% loci_missingness()
#'
#' # For missingness per locus per population
#' example_gt %>%
#'   group_by(population) %>%
#'   loci_missingness()
#' # alternatively, return a list of populations with their missingness
#' example_gt %>%
#'   group_by(population) %>%
#'   loci_missingness(type = "list")
#' # or a matrix with populations in columns and loci in rows
#' example_gt %>%
#'   group_by(population) %>%
#'   loci_missingness(type = "matrix")
#' # or within reframe (not recommended, as it much less efficient
#' # than using it directly as shown above)
#' example_gt %>%
#'   group_by(population) %>%
#'   reframe(missing = loci_missingness(genotypes))
#'
loci_missingness <- function(.x, .col = "genotypes", as_counts = FALSE,
                             n_cores = bigstatsr::nb_cores(),
                             block_size, type, ...) {
  UseMethod("loci_missingness", .x)
}

#' @export
#' @rdname loci_missingness
loci_missingness.tbl_df <- function(
    .x,
    .col = "genotypes",
    as_counts = FALSE,
    n_cores = n_cores,
    # the bigapply that splits in blocks is not
    # multithreaded, as we use the multiple
    # threads for openMP,
    block_size = bigstatsr::block_size(nrow(.x), 1), # nolint
    ...) {
  stopifnot_gen_tibble(.x)
  .col <- rlang::enquo(.col) %>%
    rlang::quo_get_expr() %>%
    rlang::as_string()
  # confirm that .col is "genotypes"
  if (.col != "genotypes") {
    stop("loci_missingness only works with the genotypes column")
  }
  loci_missingness(
    .x$genotypes,
    as_counts = as_counts,
    block_size = block_size,
    n_cores = n_cores,
    ...
  )
}


#' @export
#' @rdname loci_missingness
loci_missingness.vctrs_bigSNP <- function(
    .x,
    .col = "genotypes",
    as_counts = FALSE,
    n_cores = n_cores,
    block_size = bigstatsr::block_size(length(.x), 1), # nolint
    ...) {
  rlang::check_dots_empty()
  # get the FBM
  geno_fbm <- attr(.x, "bigsnp")$genotypes
  # rows (individuals) that we want to use
  rows_to_keep <- vctrs::vec_data(.x)
  # as long as we have more than one individual
  if (length(rows_to_keep) > 1) {
    # internal function that can be used with a big_apply
    count_na_sub <- function(BM, ind, rows_to_keep) { # nolint
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
    .col = "genotypes",
    as_counts = FALSE,
    n_cores = bigstatsr::nb_cores(),
    block_size = bigstatsr::block_size(nrow(.x), 1), # nolint
    type = c("tidy", "list", "matrix"),
    ...) {
  .col <- rlang::enquo(.col) %>%
    rlang::quo_get_expr() %>%
    rlang::as_string()
  # confirm that .col is "genotypes"
  if (.col != "genotypes") {
    stop("loci_missingness only works with the genotypes column")
  }

  # check that we only have one grouping variable
  if (length(.x %>% dplyr::group_vars()) > 1) {
    stop("loci_missingness only works with one grouping variable")
  }
  rlang::check_dots_empty()
  type <- match.arg(type)
  geno_fbm <- .gt_get_bigsnp(.x)$genotypes
  rows_to_keep <- .gt_bigsnp_rows(.x)
  count_na_sub <- function(geno_fbm, ind, rows_to_keep) {
    na_mat <- grouped_missingness_cpp( # nolint
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

  na_mat <- format_grouped_output(
    out_mat = na_mat,
    group_ids = dplyr::group_keys(.x) %>% pull(1),
    loci_names = loci_names(.x),
    type = type
  )

  return(na_mat)
}
