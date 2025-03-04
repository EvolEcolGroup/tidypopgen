#' Augment the loci table with information from a analysis object
#'
#' `augment_loci` add columns to the loci table of a `gen_tibble` related to
#' information from a given analysis.
#' @param x  An object returned by one of the `gt_` functions (e.g. [gt_pca()]).
#' @param data the `gen_tibble` used to run the PCA.
#' @param ... Additional parameters passed to the individual methods.
#' @return A [gen_tibble] with additional columns added to the loci tibble
#'   (accessible with [show_loci()]. If `data` is missing, a tibble of the
#'   information, with a column `.rownames` giving the loci names.
#' @export
augment_loci <- function(x, data, ...) {
  UseMethod("augment_loci", x)
}
