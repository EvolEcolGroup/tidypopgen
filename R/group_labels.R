#' Labels of groups, as tidypopgen names them
#'
#' This function is designed for developers of packages that add methods to
#' the generics of tidypopgen for grouped data (e.g. for genotype
#' likelihoods), so that their outputs label the groups as the methods of
#' tidypopgen do. A group is labelled by the value of its grouping variable
#' or, with several grouping variables, by their values pasted with "." (as
#' the names of the pairs of groups, e.g. `fst_pop1.pop2`).
#'
#' @param keys a data frame with one row per group and one column per
#'   grouping variable, as returned by [dplyr::group_keys()].
#' @returns a character vector with the label of each group.
#' @seealso [group_columns()] for the columns that identify the groups in a
#'   tidy output.
#' @keywords internal
#' @export
#' @examples
#' keys <- data.frame(population = c("a", "a"), sex = c("F", "M"))
#' group_labels(keys)
group_labels <- function(keys) {
  .group_labels(keys)
}
