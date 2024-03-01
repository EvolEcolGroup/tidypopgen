#' The `select_if` verb for `loci`
#'
#' An equivalent to [dplyr::select_if()] that works on the `genotype` column
#' of a `gen_tibble`. Tidy evaluation should work as expected.
#' @param .data a `gen_tibble`
#' @param .sel_arg one unquoted expression
#' @param .swap_arg not used yet
#' @param .swap_if_arg not used yet
#' @returns a list of `SNPbin` object that have been subsetted. Variable names
#' can be used as if they were positions in the data frame, so expressions
#' like x:y can be used to select a range of variables.
#' @export
#'
select_loci <-function(.data, .sel_arg, .swap_arg = NULL, .swap_if_arg = NULL){
  if (all(!is.null(.swap_arg),!is.null(.swap_if_arg))){
    stop("only one of 'swap' or 'swap_if' can be used at any one time")
  }
  # defuse the selection criteria
  sel_defused <- rlang::enquo(.sel_arg)
  # create a named vector of all loci
  loci_names_vec <- stats::setNames(seq_len(nrow(show_loci(.data))),
                             nm = show_loci_names(.data))
  # and now evaluate the selection criteria, giving it the loci names to work with
  loci_sel <- tidyselect::eval_select(expr=sel_defused, data = loci_names_vec)
  #extract the loci table
  loci_info <- attr(.data$genotypes,"loci")
  # subset the genotypes
  .data$genotypes <- lapply(.data$genotypes, .SNPbin_subset, loci_sel)
  # reintroduce the subsetted table
  attr(.data$genotypes,"loci") <- loci_info[loci_sel,]
  .data
}
