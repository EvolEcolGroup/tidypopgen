#' Count the number of loci in a `gen_tibble`
#'
#' Count the number of loci in `gen_tibble` (or directly from its `genotype`
#' column).
#' @param .x a [`gen_tibble`], or a vector of class `vctrs_bigSNP` (usually the
#'   `genotype` column of a [`gen_tibble`] object).
#' @param ... currently unused.
#' @returns the number of loci
#' @rdname count_loci
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' example_gt %>% count_loci()
count_loci <- function(.x, ...) {
  UseMethod("count_loci", .x)
}

#' @export
#' @rdname count_loci
count_loci.tbl_df <- function(.x, ...) {
  stopifnot_gen_tibble(.x)
  # extract the column and hand it over to its method
  count_loci(.x$genotypes)
}

#' @export
#' @rdname count_loci
count_loci.vctrs_bigSNP <- function(.x, ...) {
  rlang::check_dots_empty()
  nrow(attr(.x, "loci"))
}
