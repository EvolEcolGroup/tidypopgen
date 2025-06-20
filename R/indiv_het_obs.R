#' Estimate individual observed heterozygosity
#'
#' Estimate observed heterozygosity (H_obs) for each individual (i.e. the
#' frequency of loci that are heterozygous in an individual).
#'
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotype` column of
#'   a [`gen_tibble`] object), or a [`gen_tibble`].
#' @param as_counts logical, if `TRUE`, return a matrix with two columns:
#'  the number of heterozygotes and the number of missing values for each
#'  individual. These quantities can be useful to compute more complex
#'  quantities.
#' @param ... currently unused.
#' @returns either:
#'  - a vector of heterozygosities, one per individuals in the [`gen_tibble`]
#'  - a matrix with two columns, where the first is the number of heterozygous
#'    loci for each individual and the second is the number of missing values
#'    for each individual
#' @rdname indiv_het_obs
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' example_gt %>% indiv_het_obs()
#'
#' # For observed heterozygosity as counts:
#' example_gt %>% indiv_het_obs(as_counts = TRUE)
#'
indiv_het_obs <- function(.x, as_counts = FALSE, ...) {
  UseMethod("indiv_het_obs", .x)
}

#' @export
#' @rdname indiv_het_obs
indiv_het_obs.tbl_df <- function(.x, as_counts = FALSE, ...) {
  stopifnot_gen_tibble(.x)
  # extract the column and hand it over to its method
  indiv_het_obs(.x$genotypes, as_counts = as_counts, ...)
}

#' @export
#' @rdname indiv_het_obs
indiv_het_obs.vctrs_bigSNP <- function(.x, as_counts = FALSE, ...) {
  rlang::check_dots_empty()
  stopifnot_diploid(.x)
  # get the FBM
  X <- attr(.x, "bigsnp")$genotypes # nolint
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
  this_col_1_na <- bigstatsr::big_apply(
    X,
    a.FUN = count_1_na,
    ind = attr(.x, "loci")$big_index,
    a.combine = "plus",
    rows_to_keep = rows_to_keep
  )
  if (!as_counts) {
    return(this_col_1_na[1, ] / (ncol(X) - this_col_1_na[2, ]))
  } else {
    this_col_1_na <- t(this_col_1_na)
    colnames(this_col_1_na) <- c("het_n", "na_n")
    return(this_col_1_na)
  }
}

# #' @export #nolint start
# #' @rdname indiv_het_obs
# indiv_het_obs.grouped_df <- function(.x, ...){
#   .x %>% mutate(indiv_het_obs = indiv_het_obs(.data$genotypes)) %>% summarise(het_obs = mean(indiv_het_obs))
# } #nolint end
