#' Get the chromosomes of loci in a `gen_tibble`
#'
#' Extract the loci chromosomes from a `gen_tibble` (or directly from its
#' `genotype` column).
#' @param .x a [`gen_tibble`], or a vector of class `vctrs_bigSNP` (usually the
#'   `genotype` column of a [`gen_tibble`] object).
#' @param .col the column to be used when a tibble (or grouped tibble is passed
#' directly to the function). This defaults to "genotypes" and can only take
#' that value. There is no need for the user to set it, but it is included to
#' resolve certain tidyselect operations.
#' @param ... currently unused.
#' @returns a character vector of chromosomes
#' @rdname loci_chromosomes
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#' example_gt %>% loci_chromosomes()
loci_chromosomes <- function(.x, .col = "genotypes", ...) {
  UseMethod("loci_chromosomes", .x)
}

#' @export
#' @rdname loci_chromosomes
loci_chromosomes.tbl_df <- function(.x, .col = "genotypes", ...) {
  stopifnot_gen_tibble(.x)
  .col <- rlang::enquo(.col) %>%
    rlang::quo_get_expr() %>%
    rlang::as_string()
  # confirm that .col is "genotypes"
  if (.col != "genotypes") {
    stop("loci_chromosomes only works with the genotypes column")
  }
  # extract the column and hand it over to its method
  loci_chromosomes(.x$genotypes)
}

#' @export
#' @rdname loci_chromosomes
loci_chromosomes.vctrs_bigSNP <- function(.x, .col = "genotypes", ...) {
  rlang::check_dots_empty()
  attr(.x, "loci")$chromosome
}
