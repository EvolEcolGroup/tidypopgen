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
#' @rdname ind_missingness
#' @export
ind_missingness <- function(.x, as_counts = FALSE, ...) {
  UseMethod("ind_missingness", .x)
}

#' @export
#' @rdname ind_missingness
ind_missingness.tbl_df <- function(.x, as_counts = FALSE, ...){
  stopifnot_gen_tibble(.x)
  # extract the column and hand it over to its method
  ind_missingness(.x$genotypes, ...)
}

#' @export
#' @rdname ind_missingness
ind_missingness.vctrs_bigSNP <- function(.x, as_counts = FALSE, ...){
  warning("this function is not finished!!!!")
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
  this_col_na <- bigstatsr::big_apply(X, a.FUN = col_NA,
                      ind=attr(.x,"loci")$big_index,
                       a.combine = 'plus', rows_to_keep=rows_to_keep)
  if (!as_counts){
    this_col_na <- this_col_na/length(show_loci_names(.x))
  }
  this_col_na
}

#' @export
#' @rdname ind_missingness
ind_missingness.grouped_df <- function(.x, as_counts = FALSE, ...){
  .x %>% mutate(ind_missingness = ind_missingness(.data$genotypes)) %>% summarise(het_obs = mean(ind_missingness))
}
