#' The `select` verb for `loci`
#'
#' An equivalent to [dplyr::select()] that works on the `genotype` column
#' of a `gen_tibble`, using the mini-grammar available for `tidyselect`. The
#' `select`-like evaluation only has access to the names of the loci (i.e. it
#' can select only based on names, not summary statistics of those loci; look
#' at [select_loci_if()] for that feature.
#'
#' Note that the `select_loci` verb does not modify the backing FBM files,
#' but rather it subsets the list of loci to be used stored in the `gen_tibble`.
#' @param .data a `gen_tibble`
#' @param .sel_arg one unquoted expression, using the mini-grammar of
#'  [dplyr::select()] to select loci. Variable names
#' can be used as if they were positions in the data frame, so expressions
#' like x:y can be used to select a range of variables.
#' @returns a `gen_tibble` with a subset of the loci.
#' @seealso [dplyr::select()]
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # Select loci by name
#' example_gt_subset <- example_gt %>%
#'   select_loci(all_of(c("rs1", "rs2", "rs3")))
#' show_loci(example_gt_subset)
#'
#' # Select loci by index
#' example_gt_subset <- example_gt %>% select_loci(all_of(c(4, 2, 1)))
#' show_loci(example_gt_subset)
#'
select_loci <- function(.data, .sel_arg) {
  # defuse the selection criteria
  sel_defused <- rlang::enquo(.sel_arg)

  # create a named vector of all loci
  loci_names_vec <- stats::setNames(
    seq_len(nrow(show_loci(.data))),
    nm = loci_names(.data)
  )
  # now evaluate the selection criteria, giving it the loci names to work with
  loci_sel <- tidyselect::eval_select(expr = sel_defused, data = loci_names_vec)

  # modify the table
  attr(.data$genotypes, "loci") <- attr(.data$genotypes, "loci")[loci_sel, ]
  # prioritise "gen_tbl" class over "sf"
  obj_class <- class(.data)
  if ("sf" %in% obj_class) {
    obj_class <-
      c("gen_tbl", "sf", obj_class[!obj_class %in% c("gen_tbl", "sf")])
    class(.data) <- obj_class
  }
  .data
}
