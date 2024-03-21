#' Estimate individual missingness
#'
#' Estimate missingnes for each individual (i.e. the frequency of
#' missing genotypes in an individual).
#'
#' @param .x a list of [`adegenet::SNPbin`] objects (usually the `genotype` column of
#' a [`gen_tibble`] object),
#' or a [`gen_tibble`].
#' @param as_counts booelean defining whether the count of NAs (rather than the rate)
#' should be returned. It defaults to FALSE (i.e. rates are returned by default).
#' @param ... currently unused.
#' @returns a vector of heterozygosities, one per individuals in the [`gen_tibble`]
#' @rdname indiv_missingness
#' @export
indiv_missingness <- function(.x, as_counts = FALSE, ...) {
  UseMethod("indiv_missingness", .x)
}

#' @export
#' @rdname indiv_missingness
indiv_missingness.tbl_df <- function(.x, as_counts = FALSE, ...){
  stopifnot_gen_tibble(.x)
  # extract the column and hand it over to its method
  indiv_missingness(.x$genotypes, as_counts = as_counts, ...)
}

#' @export
#' @rdname indiv_missingness
indiv_missingness.vctrs_bigSNP <- function(.x, as_counts = FALSE, ...){
  rlang::check_dots_empty()
  # get the FBM
  X <- attr(.x,"bigsnp")$genotypes
  # rows (individuals) that we want to use
  rows_to_keep <- vctrs::vec_data(.x)
  # col means for submatrix (all rows, only some columns)

  # returns a matrix of 2 rows (count_1,count_na) and n_individuals columns
  col_NA <- function(X, ind, rows_to_keep) {
    count_na <- function(a){sum(is.na(a))}
    res <- apply(X[rows_to_keep,ind],1,count_na)
  }

  # count nas in one go
  this_row_na <- bigstatsr::big_apply(X, a.FUN = col_NA,
                      ind=attr(.x,"loci")$big_index,
                       a.combine = 'plus', rows_to_keep=rows_to_keep)
  if (!as_counts){
    this_row_na <- this_row_na/length(show_loci_names(.x))
  }
  this_row_na
}

#' @export
#' @rdname indiv_missingness
indiv_missingness.grouped_df <- function(.x, as_counts = FALSE, ...){
  .x %>% mutate(indiv_missingness = indiv_missingness(.data$genotypes, as_counts = as_counts)) %>%
    summarise(het_obs = mean(indiv_missingness))
}
