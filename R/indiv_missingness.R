#' Estimate individual missingness
#'
#' Estimate missingnes for each individual (i.e. the frequency of
#' missing genotypes in an individual).
#'
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotype` column of
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
  # for polyploids, this can generate a very large matrix
  # it would be better to just write a C function that counts na
  row_counts <- bigstatsr::big_counts(X, ind.col=attr(.x,"loci")$big_index,
                        ind.row = vctrs::vec_data(.x),
                        byrow = TRUE)
  row_na <- row_counts[nrow(row_counts),]
  if (!as_counts){
    row_na <- row_na/count_loci(.x)
  }
  row_na
}

#' @export
#' @rdname indiv_missingness
indiv_missingness.grouped_df <- function(.x, as_counts = FALSE, ...){
  .x %>% mutate(indiv_missingness = indiv_missingness(.data$genotypes, as_counts = as_counts)) %>%
    summarise(het_obs = mean(indiv_missingness))
}
