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
#' @param .swap_arg one unquoted expression, using the mini-grammar of
#'  [dplyr::select()] determining at which loci the reference and alterante
#'  allele will be swapped.
#' @param .swap_if_arg a logical vector of length equal to the number of loci,
#' or an expression that will tidy evaluate to such a vector.
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
  # defuse the swap criteria
  swap_defused <- rlang::enquo(.swap_arg)
  swap_if_defused <- rlang::enquo(.swap_if_arg)

  if (all(!rlang::quo_is_null(swap_defused),!rlang::quo_is_null(swap_if_defused))){
    stop("only one of 'swap' or 'swap_if' can be used at any one time")
  }
  # create a named vector of all loci
  loci_names_vec <- stats::setNames(seq_len(nrow(show_loci(.data))),
                             nm = show_loci_names(.data))
  # and now evaluate the selection criteria, giving it the loci names to work with
  loci_sel <- tidyselect::eval_select(expr=sel_defused, data = loci_names_vec)
  #extract the loci table
  loci_info <- attr(.data$genotypes,"loci")

  # subset the genotypes
  if (all(rlang::quo_is_null(swap_defused),rlang::quo_is_null(swap_if_defused))){
    .data$genotypes <- lapply(.data$genotypes, .SNPbin_subset, loci_sel)
  } else{ # subset and swap genotypes (slower)

    # if we also want to swap based on a select argument
    if (!rlang::quo_is_null(swap_defused)){
      loci_swap <- tidyselect::eval_select(expr=swap_defused, data = loci_names_vec)
      # turn them into a boolean
      #browser()
    } else { # based on an if (boolean) argument
      loci_swap <- rlang::eval_tidy(swap_if_defused,data=.data)
    }
    # swap alleles in the snpbin_list
    .data$genotypes <- lapply(.data$genotypes, select_swap_snpbin,
                              sel_indices = loci_sel,
                              to_swap = loci_swap)
    # swap alleles in table
    old_ref <- loci_info[loci_swap,"allele_ref"]
    loci_info[loci_swap,"allele_ref"] <- loci_info[loci_swap,"allele_alt"]
    loci_info[loci_swap,"allele_alt"] <- old_ref
  }
  # reintroduce the subsetted table
  attr(.data$genotypes,"loci") <- loci_info[loci_sel,]
  .data
}
