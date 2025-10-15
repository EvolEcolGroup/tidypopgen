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
#' @returns a vector of missingness, one per individuals in the
#'   [`gen_tibble`]
#' @rdname indiv_missingness
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' example_gt %>% indiv_missingness()
#'
#' # For missingness as counts:
#' example_gt %>% indiv_missingness(as_counts = TRUE)
#'
indiv_missingness <- function(.x, as_counts, block_size, ...) {
  UseMethod("indiv_missingness", .x)
}

#' @export
#' @rdname indiv_missingness
indiv_missingness.tbl_df <- function(
    .x,
    as_counts = FALSE,
    block_size = bigstatsr::block_size(nrow(.x), 1), # nolint
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
    block_size = bigstatsr::block_size(length(.x), 1), # nolint
    ...) {
  rlang::check_dots_empty()
  # get the FBM
  X <- attr(.x, "fbm") # nolint
  # rows (individuals) that we want to use
  rows_to_keep <- vctrs::vec_data(.x)


  # returns a matrix of 2 rows (count_1,count_na) and n_individuals columns
  count_1_na <- function(BM, ind, rows_to_keep) { # nolint
    gt_ind_hetero(
      BM = BM,
      rowInd = rows_to_keep,
      colInd = ind,
      ncores = 1 # n_cores, I have not seen any improvement with n_cores > 1
    )
  }

  # count heterozygotes and nas in one go
  row_na <- bigstatsr::big_apply(
    X,
    a.FUN = count_1_na,
    ind = attr(.x, "loci")$big_index,
    a.combine = "plus",
    block.size = block_size,
    rows_to_keep = rows_to_keep
  )[2, ] # get the second row (count_na)
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
