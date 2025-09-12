#' Show the ploidy information of a `gen_tibble`
#'
#' Extract the ploidy information from a  `gen_tibble`. NOTE that this function
#' does not return the ploidy level for each individual (that is obtained with
#' [`indiv_ploidy`]); instead, it returns an integer which is either the ploidy
#' level of all individuals (e.g. 2 indicates all individuals are diploid),
#' or a 0 to indicate mixed ploidy. The special case of -2 is used to indicate
#' the presence of pseudo-haploids (i.e. individuals with a ploidy of 2 but
#' for which we only have information for one allele; the dosages are 0 or 2
#' for these individuals).
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotype` column of
#' a [`gen_tibble`] object),
#' or a [`gen_tibble`].
#' @param ... currently unused.
#' @returns the ploidy (0 indicates mixed ploidy)
#' @rdname show_ploidy
#' @export
#' @seealso [indiv_ploidy()]
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' example_gt %>% show_ploidy()
show_ploidy <- function(.x, ...) {
  UseMethod("show_ploidy", .x)
}

#' @export
#' @rdname show_ploidy
show_ploidy.tbl_df <- function(.x, ...) {
  stopifnot_gen_tibble(.x)
  # extract the column and hand it over to its method
  show_ploidy(.x$genotypes)
}


#' @export
#' @rdname show_ploidy
show_ploidy.vctrs_bigSNP <- function(.x, ...) {
  rlang::check_dots_empty()
  attr(.x, "ploidy")
}
