#' Return individual ploidy
#'
#' Returns the ploidy for each individual. On a grouped [`gen_tibble`], the
#' ploidy of all the individuals is returned, in the order of the rows (it
#' does not depend on the groups).
#'
#' @param .x a [`gen_tibble`] (grouped or not), or a vector of class
#'   `vctrs_bigSNP` (usually the `genotype` column of a [`gen_tibble`] object)
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
  stopifnot_has_genotypes(.x)
  # extract the column and hand it over to its method
  indiv_ploidy(.x$genotypes, ...)
}

#' @export
#' @rdname indiv_ploidy
indiv_ploidy.grouped_gen_tbl <- function(.x, ...) {
  stopifnot_has_genotypes(.x)
  # ploidy does not depend on the groups: return it for all the individuals
  # (the whole column), in the order of the rows
  indiv_ploidy(.x$genotypes, ...)
}

# (no method for grouped_df: ploidy does not depend on the groups, so other
# grouped tibbles, including those of other packages built on tidypopgen, use
# the tbl_df method, which hands their genotypes column over to its method)

# a data frame that has lost the class of its tibble (e.g. a gen_tibble
# returned by pick() or reframe()): restore the class and dispatch again
#' @export
#' @noRd
indiv_ploidy.data.frame <- function(.x, ...) {
  indiv_ploidy(restore_genotypes_tbl(.x, "indiv_ploidy"), ...)
}

#' @export
#' @rdname indiv_ploidy
indiv_ploidy.vctrs_bigSNP <- function(.x, ...) {
  rlang::check_dots_empty()
  if ((show_ploidy(.x) != 0) && (show_ploidy(.x) != -2)) {
    rep(show_ploidy(.x), length(.x))
  } else {
    attr(.x, "fbm_ploidy")[vctrs::vec_data(.x)]
  }
}
