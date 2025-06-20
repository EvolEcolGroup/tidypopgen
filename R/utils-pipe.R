#' Pipe operator
#'
#' See `magrittr::pipe \%>\% for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom dplyr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#' example_gt %>% count_loci()
NULL
