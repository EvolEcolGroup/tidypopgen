#' Find transitions
#'
#' Use the loci table to define which loci are transitions
#'
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotype` column of
#' a [`gen_tibble`] object),
#' or a [`gen_tibble`].
#' @param .col the column to be used when a tibble (or grouped tibble is passed
#' directly to the function). This defaults to "genotypes" and can only take
#' that value. There is no need for the user to set it, but it is included to
#' resolve certain tidyselect operations.
#' @param ... other arguments passed to specific methods.
#' @returns a logical vector defining which loci are transitions
#' @rdname loci_transitions
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#' example_gt %>% loci_transitions()
loci_transitions <- function(.x, .col = "genotypes", ...) {
  UseMethod("loci_transitions", .x)
}

#' @export
#' @rdname loci_transitions
loci_transitions.tbl_df <- function(.x, .col = "genotypes", ...) {
  stopifnot_gen_tibble(.x)
  .col <- rlang::enquo(.col) %>%
    rlang::quo_get_expr() %>%
    rlang::as_string()
  # confirm that .col is "genotypes"
  if (.col != "genotypes") {
    stop("loci_transitions only works with the genotypes column")
  }
  check_allele_alphabet(show_loci(.x$genotypes))
  loci_transitions(.x$genotypes, ...)
}


#' @export
#' @rdname loci_transitions
loci_transitions.vctrs_bigSNP <- function(.x, .col = "genotypes", ...) {
  rlang::check_dots_empty()

  !loci_transversions(.x)
}
