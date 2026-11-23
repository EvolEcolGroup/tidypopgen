#' Columns that identify groups in a tidy output, as tidypopgen names them
#'
#' This function is designed for developers of packages that add methods to
#' the generics of tidypopgen for grouped data (e.g. for genotype
#' likelihoods), so that their tidy outputs identify the groups as the
#' methods of tidypopgen do: a single `group` column when there is one
#' grouping variable, and one column per grouping variable when there are
#' several (as [dplyr::summarise()]), keeping their types (e.g. a factor keeps
#' its levels). For outputs on pairs of groups, each group of the pair has one
#' column per grouping variable, named with a suffix (e.g. `population_1` and
#' `population_2`, with `suffix = "_1"` and `"_2"`), also when there is a
#' single grouping variable.
#'
#' @param keys a data frame with one row per group and one column per
#'   grouping variable, as returned by [dplyr::group_keys()].
#' @param idx the index of the group of each row of the output.
#' @param suffix `NULL` for outputs with a row per group, or a suffix added to
#'   the names of the grouping variables, for the columns of one group of a
#'   pair.
#' @param other_names the names of the other columns of the output: the group
#'   columns must not have any of these names (an error asks to rename the
#'   grouping variables).
#' @returns a tibble with one row per element of `idx`.
#' @seealso [group_labels()] for the labels of the groups.
#' @keywords internal
#' @export
#' @examples
#' keys <- data.frame(population = c("a", "a", "b"), sex = c("F", "M", "F"))
#' group_columns(keys, c(1, 3))
#' group_columns(keys[1], c(1, 3))
#' group_columns(keys, c(1, 3), suffix = "_1")
group_columns <- function(keys, idx, suffix = NULL, other_names = character()) {
  keys <- tibble::as_tibble(keys)
  if (is.null(suffix)) {
    cols <- .group_columns(list(vars = names(keys), keys = keys), idx)
  } else {
    cols <- keys[idx, , drop = FALSE]
    names(cols) <- paste0(names(keys), suffix)
  }
  .check_group_col_names(names(cols), other_names)
  cols
}
