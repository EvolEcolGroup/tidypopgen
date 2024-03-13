#' Estimate individual observed heterozygosity
#'
#' Estimate observed heterozygosity (H_obs) for each individual (i.e. the frequency of
#' loci that are heterozygous in an individual).
#'
#' @param .x a list of [`adegenet::SNPbin`] objects (usually the `genotype` column of
#' a [`gen_tibble`] object),
#' or a [`gen_tibble`].
#' @param ... currently unused.
#' @returns a vector of heterozygosities, one per individuals in the [`gen_tibble`]
#' @rdname ind_het_obs
#' @export
ind_het_obs <- function(.x, ...) {
  UseMethod("ind_het_obs", .x)
}

#' @export
#' @rdname ind_het_obs
ind_het_obs.tbl_df <- function(.x, ...){
  stopifnot_gen_tibble(.x)
  # extract the column and hand it over to its method
  ind_het_obs(.x$genotypes, ...)
}

#' @export
#' @rdname ind_het_obs
ind_het_obs.vctrs_bigSNP <- function(.x, ...){
  rlang::check_dots_empty()
  X <- attr(.x,"bigsnp")$genotypes


  rowsums <- bigstatsr::big_apply(X, a.FUN = function(X, ind) row_count_1(X[, ind]),
                      ind=attr(.x,"loci")$big_index,
                       a.combine = 'plus')
  rowNA <- bigstatsr::big_apply(X, a.FUN = function(X, ind) row_count_NA(X[, ind]),
                                ind=attr(.x,"loci")$big_index,
                     a.combine = 'plus')
  rowsums/(ncol(X)-rowNA)
}


row_count_NA <- function(X){
  count_na <- function(a){sum(is.na(a))}
  apply(X,1,count_na)
}

row_count_1 <- function(X){
  count_1 <- function(a){sum(a==1, na.rm = TRUE)}
  apply(X,1,count_1)
}
