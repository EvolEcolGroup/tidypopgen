#' Get the names of loci in a `gen_tibble`
#'
#' Extract the loci names from a  `gen_tibble` (or directly from its `genotype`
#' column).
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotype` column of
#' a [`gen_tibble`] object),
#' or a [`gen_tibble`].
#' @param ... currently unused.
#' @returns a character vector of names
#' @rdname loci_names
#' @export
#' @examples
#' example_gt <- example_gt()
#' example_gt %>% loci_names()
loci_names <- function(.x, ...) {
  UseMethod("loci_names", .x)
}

#' @export
#' @rdname loci_names
loci_names.tbl_df <- function(.x, ...) {
  stopifnot_gen_tibble(.x)
  # extract the column and hand it over to its method
  loci_names(.x$genotypes)
}

#' @export
#' @rdname loci_names
loci_names.vctrs_bigSNP <- function(.x, ...) {
  rlang::check_dots_empty()
  attr(.x, "loci")$name
}
