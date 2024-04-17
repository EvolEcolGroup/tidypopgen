#' The `select_if` verb for `loci`
#'
#' An equivalent to [dplyr::select_if()] that works on the `genotype` column
#' of a `gen_tibble`. This function has access to the genotypes (and thus can
#' work on summary statistics to select), but not the names of the loci (look
#' at [select_loci()] for that feature.
#'
#' #' Note that the `select_loci_if` verb does not modify the backing FBM files,
#' but rather it subsets the list of loci to be used stored in the `gen_tibble`.
#' @param .data a `gen_tibble`
#' @param .sel_logical a logical vector of length equal to the number of loci,
#' or an expression that will tidy evaluate to such a vector.
#' @export
#'
select_loci_if <-function(.data, .sel_logical){
  # defuse the boolean argument
  sel_defused <- rlang::enquo(.sel_logical)

  # and now evaluate it, allowing it to see the data
  loci_sel <- rlang::eval_tidy(sel_defused,data=.data)
  if (!inherits(loci_sel,"logical")){
    stop(".sel_logical should be a logical (boolean) vector")
  }
  if (length(loci_sel) != ncol(show_genotypes(.data$genotypes))){
    stop(".sel_logical should be the same length as the number of loci")
  }

  # subset the loci table
  attr(.data$genotypes,"loci") <- attr(.data$genotypes,"loci")[loci_sel,]
  .data
}

