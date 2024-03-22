#' Count the number of loci in a `gen_tibble`
#'
#' Count the number of loci in `gen_tibble` (or directly from its `genotype`
#' column).
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotype` column of
#' a [`gen_tibble`] object),
#' or a [`gen_tibble`].
#' @param ... currently unused.
#' @returns a character vector of names
#' @rdname count_loci
#' @export
count_loci <- function(.x, ...) {
  UseMethod("count_loci", .x)
}

#' @export
#' @rdname count_loci
count_loci.tbl_df <- function(.x, ...){
  stopifnot_gen_tibble(.x)
  # extract the column and hand it over to its method
  count_loci(.x$genotypes)
}

#' @export
#' @rdname count_loci
count_loci.vctrs_bigSNP <- function(.x, ...){
  rlang::check_dots_empty()
  nrow(attr(.x,"loci"))
}
