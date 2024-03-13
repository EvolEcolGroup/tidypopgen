#' The `select` verb for `loci`
#'
#' An equivalent to [dplyr::select()] that works on the `genotype` column
#' of a `gen_tibble`, using the mini-grammar available for `tidyselect`. The
#' `select`-like evaluation only has access to the names of the loci (i.e. it
#' can select only based on names, not summary statistics of those loci; look
#' at [select_loci_if()] for that feature. There is also the possibility
#' of swapping reference and alternate alleles whilst selecting loci, using
#' either the argument `.swap_arg` (which has `select`-like mini-grammar on names)
#' or `.swap_arg_if` (which has access to the genotypes like `select_loci_if`)
#' @param .data a `gen_tibble`
#' @param .sel_arg one unquoted expression, using the mini-grammar of
#'  [dplyr::select()] to select loci
#' @returns a list of `SNPbin` object that have been sub-setted to the selected
#' loci.
#'
#' Variable names
#' can be used as if they were positions in the data frame, so expressions
#' like x:y can be used to select a range of variables.
#' @export
#'
select_loci <-function(.data, .sel_arg, .swap_arg = NULL, .swap_if_arg = NULL){
  # defuse the selection criteria
  sel_defused <- rlang::enquo(.sel_arg)

  # create a named vector of all loci
  loci_names_vec <- stats::setNames(seq_len(nrow(show_loci(.data))),
                             nm = show_loci_names(.data))
  # and now evaluate the selection criteria, giving it the loci names to work with
  loci_sel <- tidyselect::eval_select(expr=sel_defused, data = loci_names_vec)

  # modify the table
  attr(.data$genotypes,"loci") <- attr(.data$genotypes,"loci")[loci_sel,]
  .data
}
