#' Count number of non-missing genotypes at each locus
#'
#' Count number of non-missing genotypes at each locus. Note that this function returns a count
#' of genotypes, not a count of alleles (which is twice the number of genotypes).
#'
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotypes` column of
#' a [`gen_tibble`] object),
#' or a [`gen_tibble`].
#' @param n_cores number of cores to be used, it defaults to [bigstatsr::nb_cores()]
#' @param ... other arguments passed to specific methods.
#' @returns a vector of frequencies, one per locus
#' @rdname loci_n
#' @export
loci_n <- function(.x, ...) {
  UseMethod("loci_n", .x)
}

#' @export
#' @rdname loci_n
loci_n.tbl_df <- function(.x, ...) {
  # TODO this is a hack to deal with the class being dropped when going through group_map
  stopifnot_gen_tibble(.x)
  loci_n(.x$genotypes, ...)
}


#' @export
#' @rdname loci_n
loci_n.vctrs_bigSNP <- function(.x, ...) {
  rlang::check_dots_empty()

  tot_n <- length(.x)
  #  stopifnot_diploid(.x)
  loc_n <- tot_n - .x %>% loci_missingness(as_counts = TRUE)
  loc_n
}

#' @export
#' @rdname loci_n
loci_n.grouped_df <- function(.x, n_cores = bigstatsr::nb_cores(), ...) {
  rlang::check_dots_empty()

  # return a list to mimic a group_map
  lapply(seq_len(ncol(na_mat)), function(i) na_mat[, i])
}
