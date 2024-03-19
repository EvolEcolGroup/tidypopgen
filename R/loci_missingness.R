#' Estimate missingness at each locus
#'
#' Estimate the rate of missingness at each locus.
#'
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotypes` column of
#' a [`gen_tibble`] object),
#' or a [`gen_tibble`].
#' @param as_counts booelean defining whether the count of NAs (rather than the rate)
#' should be returned. It defaults to FALSE (i.e. rates are returned by default).
#' @param ... other arguments passed to specific methods.
#' @returns a vector of frequencies, one per locus
#' @rdname loci_missingness
#' @export
loci_missingness <- function(.x, as_counts = FALSE, ...) {
  UseMethod("loci_missingness", .x)
}

# We should write a cpp counts function. We can't use the snp_* family of functions
# as they ignore NAs

#' @export
#' @rdname loci_missingness
loci_missingness.tbl_df <- function(.x, as_counts = FALSE, ...) {
  #TODO this is a hack to deal with the class being dropped when going through group_map
  stopifnot_gen_tibble(.x)
  loci_missingness(.x$genotypes, as_counts = as_counts, ...)
}


#' @export
#' @rdname loci_missingness
loci_missingness.vctrs_bigSNP <- function(.x, as_counts = FALSE, ...) {
  rlang::check_dots_empty()
  # get the FBM
  geno_fbm <- attr(.x,"bigsnp")$genotypes
  # rows (individuals) that we want to use
  rows_to_keep <- vctrs::vec_data(.x)
  # as long as we have more than one individual
  if (length(rows_to_keep)>1){
    # # col means for submatrix (all rows, only some columns)
    # colMeans_sub <- function(X, ind, rows_to_keep) {
    #   count_na <- function(x){sum(is.na(x))}
    #   apply(X[rows_to_keep, ind], 2, count_na)
    # }
    # n_na <- bigstatsr::big_apply(geno_fbm, a.FUN = colMeans_sub,
    #                              rows_to_keep = rows_to_keep,
    #                              ind=attr(.x,"loci")$big_index,
    #                              a.combine = 'c')
    n_na <- bigstatsr::big_counts(geno_fbm, ind.row = rows_to_keep,
                          ind.col = attr(.x,"loci")$big_index)[4,]
    if (!as_counts){
      n_na <- n_na/length(rows_to_keep)
    }
  } else { # if we have a single individual
    n_na <-geno_fbm[rows_to_keep,attr(.x,"loci")$big_index]
  }
  n_na
}

#' @export
#' @rdname loci_missingness
loci_missingness.grouped_df <- function(.x, as_counts = FALSE, ...) {
  # TODO this is seriously inefficient, we need to cast it into a big_apply problem
  # of maybe it isn't that bad...
  group_map(.x, .f=~loci_missingness(.x, as_counts=as_counts))
}

