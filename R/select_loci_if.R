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
#' or an expression that will tidy evaluate to such a vector. Only loci
#' for which .sel_logical is TRUE will be selected; NA will be treated as FALSE.
#' @returns a subset of the list of loci in the `gen_tibble`
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # Select loci by chromosome
#' example_gt_subset <- example_gt %>%
#'   select_loci_if(loci_chromosomes(genotypes) == "chr2")
#' show_loci(example_gt_subset)
#'
#' # Select loci by a summary statistic
#' example_gt_subset <- example_gt %>%
#'   select_loci_if(loci_maf(genotypes) > 0.2)
#' show_loci(example_gt_subset)
#'
select_loci_if <- function(.data, .sel_logical) {
  # defuse the boolean argument
  sel_defused <- rlang::enquo(.sel_logical)

  # and now evaluate it, allowing it to see the data
  loci_sel <- rlang::eval_tidy(sel_defused, data = .data)
  if (!inherits(loci_sel, "logical")) {
    stop(".sel_logical should be a logical boolean vector")
  }
  if (length(loci_sel) != ncol(show_genotypes(.data$genotypes))) {
    stop(".sel_logical should be the same length as the number of loci")
  }
  # set to FALSE any missing value
  loci_sel[is.na(loci_sel)] <- FALSE

  # subset the loci table
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
