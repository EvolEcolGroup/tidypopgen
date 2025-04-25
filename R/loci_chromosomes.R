#' Get the chromosomes of loci in a `gen_tibble`
#'
#' Extract the loci chromosomes from a `gen_tibble` (or directly from its
#' `genotype` column).
#' @param .x a [`gen_tibble`], or a vector of class `vctrs_bigSNP` (usually the
#'   `genotype` column of a [`gen_tibble`] object).
#' @param ... currently unused.
#' @returns a character vector of chromosomes
#' @rdname loci_chromosomes
#' @export
#' @examples
#' example_gt <- example_gen_tibble()
#' example_gt %>% loci_chromosomes()
loci_chromosomes <- function(.x, ...) {
  UseMethod("loci_chromosomes", .x)
}

#' @export
#' @rdname loci_chromosomes
loci_chromosomes.tbl_df <- function(.x, ...) {
  stopifnot_gen_tibble(.x)
  # extract the column and hand it over to its method
  loci_chromosomes(.x$genotypes)
}

#' @export
#' @rdname loci_chromosomes
loci_chromosomes.vctrs_bigSNP <- function(.x, ...) {
  rlang::check_dots_empty()
  attr(.x, "loci")$chromosome
}
