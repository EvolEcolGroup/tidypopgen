#' Restore the class of a tibble from its genotypes column
#'
#' Some `dplyr` operations return tibbles that have lost the class of a
#' `gen_tibble`, such as the data frames returned by `pick()` or `across()`
#' within `summarise()`, `reframe()` or `mutate()`, the result of `reframe()`,
#' or `as_tibble()`. Their `genotypes` column, however, keeps its class.
#' The generics of `tidypopgen` therefore have fallback methods for data
#' frames (and grouped data frames), which use `restore_genotypes_tbl()` to
#' give back the tibble its class, as given by `genotypes_tbl_class()` for
#' its `genotypes` column, and then dispatch again.
#'
#' Packages that build their own tibbles on `tidypopgen` (e.g. for genotype
#' likelihoods) add a method to `genotypes_tbl_class()` for the
#' class of their `genotypes` column, so that their own methods of the
#' `tidypopgen` generics are used after the class is restored.
#'
#' @param genotypes the `genotypes` column of a tibble
#' @param ... currently unused
#' @returns `genotypes_tbl_class()` returns the class of the tibble that holds
#'   the `genotypes` column (`"gen_tbl"` for a `vctrs_bigSNP` column), or
#'   `NULL` for a column of an unknown class.
#' @keywords internal
#' @export
genotypes_tbl_class <- function(genotypes, ...) {
  UseMethod("genotypes_tbl_class", genotypes)
}

#' @export
#' @rdname genotypes_tbl_class
genotypes_tbl_class.default <- function(genotypes, ...) {
  NULL
}

#' @export
#' @rdname genotypes_tbl_class
genotypes_tbl_class.vctrs_bigSNP <- function(genotypes, ...) {
  "gen_tbl"
}

#' @param .x a data frame or tibble with a `genotypes` column
#' @param generic the name of the generic for which the class is restored,
#'   used in the error message if the restored class has no method
#' @returns `restore_genotypes_tbl()` returns `.x` as a tibble with the class
#'   given by `genotypes_tbl_class()`, and its grouped version (e.g.
#'   `"grouped_gen_tbl"`) if `.x` is grouped. A `rowwise_df` is ungrouped. It
#'   stops with an error if `.x` already has that class, as the generic would
#'   otherwise dispatch to the fallback method again.
#' @rdname genotypes_tbl_class
#' @export
restore_genotypes_tbl <- function(.x, generic = "this function") {
  if (!"genotypes" %in% names(.x)) {
    stop("not a gen_tibble, 'genotypes' column is missing")
  }
  tbl_class <- genotypes_tbl_class(.x$genotypes)
  if (is.null(tbl_class)) {
    stop("not a gen_tibble, the genotypes column is not of class vctrs_bigSNP")
  }
  if (inherits(.x, "rowwise_df")) {
    .x <- dplyr::ungroup(.x)
  }
  if (!inherits(.x, "tbl_df")) {
    .x <- tibble::as_tibble(.x)
  }
  if (inherits(.x, "grouped_df")) {
    new_class <- c(paste0("grouped_", tbl_class), "grouped_df", tbl_class)
  } else {
    new_class <- tbl_class
  }
  if (inherits(.x, new_class[1])) {
    stop(generic, "() has no method for objects of class '", new_class[1], "'")
  }
  # this only copies the list of columns, not the columns themselves
  class(.x) <- unique(c(new_class, class(.x)))
  .x
}
