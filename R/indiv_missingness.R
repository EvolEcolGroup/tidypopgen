#' Estimate individual missingness
#'
#' Estimate missingness for each individual (i.e. the frequency of missing
#' genotypes in an individual).
#'
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotype` column of
#'   a [`gen_tibble`] object), or a [`gen_tibble`].
#' @param as_counts boolean defining whether the count of NAs (rather than the
#'   rate) should be returned. It defaults to FALSE (i.e. rates are returned by
#'   default).
#' @param block_size maximum number of loci read at once.
#' @param ... currently unused.
#' @returns a vector of heterozygosities, one per individuals in the
#'   [`gen_tibble`]
#' @rdname indiv_missingness
#' @export
indiv_missingness <- function(.x, as_counts, block_size, ...) {
  UseMethod("indiv_missingness", .x)
}

#' @export
#' @rdname indiv_missingness
indiv_missingness.tbl_df <- function(
    .x,
    as_counts = FALSE,
    block_size = bigstatsr::block_size(nrow(attr(.x, "loci")), 1), # nolint
    ...) {
  stopifnot_gen_tibble(.x)
  # extract the column and hand it over to its method
  indiv_missingness(
    .x$genotypes,
    as_counts = as_counts,
    block_size = block_size,
    ...
  )
}

#' @export
#' @rdname indiv_missingness
indiv_missingness.vctrs_bigSNP <- function(
    .x,
    as_counts = FALSE,
    block_size = bigstatsr::block_size(nrow(attr(.x, "loci")), 1), # nolint
    ...) {
  rlang::check_dots_empty()
  # get the FBM
  X <- attr(.x, "bigsnp")$genotypes # nolint
  # rows (individuals) that we want to use
  rows_to_keep <- vctrs::vec_data(.x)

  # for polyploids, this can generate a very large matrix
  # it would be better to just write a C function that counts na
  count_row_na_sub <- function(X, ind, rows_to_keep) { # nolint
    row_counts <- bigstatsr::big_counts(
      X, # nolint
      ind.col = ind,
      ind.row = rows_to_keep,
      byrow = TRUE
    )
    row_na <- row_counts[nrow(row_counts), ] # nolint
  }

  row_na <- bigstatsr::big_apply(
    X,
    a.FUN = count_row_na_sub,
    rows_to_keep = rows_to_keep,
    ind = attr(.x, "loci")$big_index,
    ncores = 1, # parallelisation is used within the function
    block.size = block_size,
    a.combine = "plus"
  )
  if (!as_counts) {
    row_na <- row_na / count_loci(.x)
  }
  row_na
}

#' #' @export
#' #' @rdname indiv_missingness
#' indiv_missingness.grouped_df <- function(.x, as_counts = FALSE, ...){
#'   .x %>% mutate(indiv_missingness = indiv_missingness(.data$genotypes,
#'   as_counts = as_counts)) %>%
#'     summarise(het_obs = mean(indiv_missingness))
#' }
