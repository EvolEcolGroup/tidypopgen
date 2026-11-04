# Internal helpers to describe the groups of a grouped gen_tibble, so that
# every function handles one or more grouping variables, and empty groups, in
# the same way.
#
# Groups are taken from dplyr::group_data(), and are thus in the order of
# dplyr::group_keys(). They include empty groups, which exist only if the user
# asked for them (with group_by(.drop = FALSE), or filter() / slice() with
# .preserve = TRUE); functions report NA (or NaN) for them.
#
# When a single label is needed for each group (e.g. names of a vector,
# column names of a matrix, a column of individual metadata), the values of
# the grouping variables are pasted with ".", as in the names of the pairs of
# groups (e.g. fst_pop1.pop2). Tidy outputs with a row per group have one
# column per grouping variable when there are several (as with
# dplyr::summarise()), and a single `group` column otherwise; either way, the
# columns keep the type of the grouping variables (e.g. a factor keeps its
# levels), so that sorting or plotting them keeps the order of the groups.

#' Information on the groups of a grouped gen_tibble
#'
#' @param .x a grouped gen_tibble
#' @returns a list with `vars` (the grouping variables), `keys` (a tibble
#'   with one row per group and one column per grouping variable, as
#'   [dplyr::group_keys()]), `n` (the number of groups, including empty
#'   ones), `ids` (the 0-based group of each individual, as needed by the cpp
#'   code), `labels` (a character label for each group) and `sizes` (the
#'   number of individuals in each group)
#' @keywords internal
#' @noRd
.group_info <- function(.x) {
  keys <- dplyr::group_keys(.x)
  list(
    vars = dplyr::group_vars(.x),
    keys = keys,
    n = nrow(keys),
    ids = dplyr::group_indices(.x) - 1L,
    labels = .group_labels(keys),
    sizes = lengths(dplyr::group_rows(.x))
  )
}

#' A character label for each group
#'
#' @param keys a tibble of group keys, as returned by [dplyr::group_keys()]
#' @returns a character vector with one label per group: the value of the
#'   grouping variable or, with several grouping variables, their values
#'   pasted with "."
#' @keywords internal
#' @noRd
.group_labels <- function(keys) {
  if (ncol(keys) == 1) {
    as.character(keys[[1]])
  } else {
    do.call(paste, c(unname(lapply(keys, as.character)), sep = "."))
  }
}

#' The group of each individual
#'
#' @param .x a grouped gen_tibble
#' @returns a vector with the group of each individual (row) of `.x`: the
#'   grouping variable itself (keeping its type) if there is only one, or
#'   the group labels (see `.group_labels()`) if there are several, as a
#'   factor with the labels as levels in the order of the groups (rather than
#'   characters, which would sort alphabetically)
#' @keywords internal
#' @noRd
.group_of_rows <- function(.x) {
  vars <- dplyr::group_vars(.x)
  if (length(vars) == 1) {
    .x[[vars]]
  } else {
    labels <- .group_labels(dplyr::group_keys(.x))
    factor(labels[dplyr::group_indices(.x)], levels = labels)
  }
}

#' Columns identifying the groups in a tidy output
#'
#' @param groups the group information, as returned by `.group_info()`
#' @param idx the index of the group of each row of the output
#' @returns a tibble with one row per element of `idx`: a single `group`
#'   column if there is one grouping variable, or one column per grouping
#'   variable if there are several, with the values of the grouping variables
#'   (keeping their type, e.g. a factor keeps its levels)
#' @keywords internal
#' @noRd
.group_columns <- function(groups, idx) {
  if (length(groups$vars) == 1) {
    tibble::new_tibble(list(group = groups$keys[[1]][idx]), nrow = length(idx))
  } else {
    groups$keys[idx, , drop = FALSE]
  }
}

#' Bind the results of each group into a tidy tibble
#'
#' The equivalent of `dplyr::bind_rows(res, .id = "group")` for a list named
#' with the group labels, with the group columns of `.group_columns()` first.
#'
#' @param res a list of data frames, one per group (in group order)
#' @param groups the group information, as returned by `.group_info()`
#' @returns a tibble
#' @keywords internal
#' @noRd
.bind_group_rows <- function(res, groups) {
  idx <- rep(seq_along(res), vapply(res, nrow, integer(1)))
  dplyr::bind_cols(
    .group_columns(groups, idx),
    dplyr::bind_rows(unname(res))
  )
}

#' Labels for pairs of groups
#'
#' @param labels the group labels, as returned by `.group_labels()`
#' @param pairwise_combn a matrix with two rows, giving the indices of the
#'   two groups of each pair (as returned by `utils::combn()`)
#' @returns a character vector with the label of each pair (e.g. pop1.pop2)
#' @keywords internal
#' @noRd
.pair_labels <- function(labels, pairwise_combn) {
  paste(labels[pairwise_combn[1, ]], labels[pairwise_combn[2, ]], sep = ".")
}
