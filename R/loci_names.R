#' Get the names of loci in a `gen_tibble`
#'
#' Extract the loci names from a  `gen_tibble` (or directly from its `genotype`
#' column).
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotype` column of
#' a [`gen_tibble`] object),
#' or a [`gen_tibble`].
#' @param .col the column to be used when a tibble (or grouped tibble is passed
#' directly to the function). This defaults to "genotypes" and can only take
#' that value. There is no need for the user to set it, but it is included to
#' resolve certain tidyselect operations.
#' @param ... currently unused.
#' @returns a character vector of names
#' @rdname loci_names
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#' example_gt %>% loci_names()
loci_names <- function(.x, .col = "genotypes", ...) {
  UseMethod("loci_names", .x)
}

#' @export
#' @rdname loci_names
loci_names.tbl_df <- function(.x, .col = "genotypes", ...) {
  stopifnot_gen_tibble(.x)
  # extract the column and hand it over to its method
  .col <- rlang::enquo(.col) %>%
    rlang::quo_get_expr() %>%
    rlang::as_string()
  # confirm that .col is "genotypes"
  if (.col != "genotypes") {
    stop("loci_names only works with the genotypes column")
  }
  loci_names(.x$genotypes)
}

#' @export
#' @rdname loci_names
loci_names.vctrs_bigSNP <- function(.x, .col = "genotypes", ...) {
  rlang::check_dots_empty()
  attr(.x, "loci")$name
}
