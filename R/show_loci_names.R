#' Show the loci names of a `gen_tibble`
#'
#' Extract the loci names from a  `gen_tibble` (or directly from its `genotype`
#' column).
#' @param .x a list of [`adegenet::SNPbin`] objects (usually the `genotype` column of
#' a [`gen_tibble`] object),
#' or a [`gen_tibble`].
#' @param ... currently unused.
#' @returns a character vector of names
#' @rdname show_loci_names
#' @export
show_loci_names <- function(.x, ...) {
  UseMethod("show_loci_names", .x)
}

#' @export
#' @rdname show_loci_names
show_loci_names.tbl_df <- function(.x, ...){
  stopifnot_gen_tibble(.x)
  # extract the column and hand it over to its method
  show_loci_names(.x$genotypes)
}

#' @export
#' @rdname show_loci_names
show_loci_names.vctrs_bigSNP <- function(.x, ...){
  rlang::check_dots_empty()
  attr(.x,"loci")$name
}
