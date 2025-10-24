#' Order the loci table of a gen_tibble
#'
#' This function reorders the loci table so that positions within a chromosome
#' are sequential. It also re-saves the genotypes into a new file backed matrix
#' with the new order, so that it can be used by functions such as
#' [loci_ld_clump()] and [gt_pca_autoSVD()]. If the loci table is already
#' ordered, the original `gen_tibble` is returned.
#'
#' @param .x a [`gen_tibble`]
#' @param use_current_table boolean, if FALSE (the default), the table will be
#'   reordered; if TRUE, then the current loci table, which might have been
#'   reordered manually, will be used, but only if the positions within each
#'   chromosome are sequential
#' @param ignore_genetic_dist boolean to ignore the genetic distance when
#'   checking. Note that, if `genetic_dist` are being ignored and they are not
#'   sorted, the function will set them to zero to avoid problems with other
#'   software.
#' @param quiet boolean to suppress information about the files
#' @param ... other arguments
#' @return A [gen_tibble]
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl") %>% select_loci(c(1, 5, 2, 6, 4, 3))
#'
#' # Loci are in the wrong order
#' show_loci(example_gt)
#'
#' # Reorder the loci, ignoring genetic distance
#' example_gt_ordered <- gt_order_loci(example_gt, ignore_genetic_dist = TRUE)
#'
#' # Loci are now in the correct order
#' show_loci(example_gt_ordered)
gt_order_loci <- function(
    .x,
    use_current_table = FALSE,
    ignore_genetic_dist = TRUE,
    quiet = FALSE,
    ...) {
  if (use_current_table) {
    new_table <- show_loci(.x)
  } else {
    new_table <- show_loci(.x) %>%
      dplyr::arrange(cast_chromosome_to_int(.data$chromosome), .data$position)
    show_loci(.x) <- new_table
  }
  # if asked to use the current table, check that it is ordered
  is_loci_table_ordered(
    .x,
    error_on_false = TRUE,
    ignore_genetic_dist = ignore_genetic_dist
  )
  gt_update_backingfile(
    .x,
    quiet = quiet,
    rm_unsorted_dist = ignore_genetic_dist,
    ...
  )
}
