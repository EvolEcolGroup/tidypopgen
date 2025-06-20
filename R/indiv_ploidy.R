#' Return individual ploidy
#'
#' Returns the ploidy for each individual.
#'
#' @param .x a [`gen_tibble`], or a vector of class `vctrs_bigSNP` (usually the
#'   `genotype` column of a [`gen_tibble`] object)
#' @param ... currently unused.
#' @returns a vector of ploidy, one per individuals in the [`gen_tibble`]
#' @rdname indiv_ploidy
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' example_gt %>% indiv_ploidy()
#'
indiv_ploidy <- function(.x, ...) {
  UseMethod("indiv_ploidy", .x)
}

#' @export
#' @rdname indiv_ploidy
indiv_ploidy.tbl_df <- function(.x, ...) {
  stopifnot_gen_tibble(.x)
  # extract the column and hand it over to its method
  indiv_ploidy(.x$genotypes, ...)
}

#' @export
#' @rdname indiv_ploidy
indiv_ploidy.vctrs_bigSNP <- function(.x, ...) {
  rlang::check_dots_empty()
  if ((show_ploidy(.x) != 0) && (show_ploidy(.x) != -2)) {
    rep(show_ploidy(.x), length(.x))
  } else {
    attr(.x, "bigsnp")$fam$ploidy[vctrs::vec_data(.x)]
  }
}

#' #' @export
#' #' @rdname indiv_ploidy
#' indiv_ploidy.grouped_df <- function(.x, ...){
#'   .x %>% mutate(indiv_ploidy = indiv_ploidy(.data$genotypes)) %>%
#'     summarise(ploidy = mean(indiv_ploidy))
#' }
