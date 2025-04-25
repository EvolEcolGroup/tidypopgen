#' Find transversions
#'
#' Use the loci table to define which loci are transversions
#'
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotype` column of
#' a [`gen_tibble`] object),
#' or a [`gen_tibble`].
#' @param ... other arguments passed to specific methods.
#' @returns a logical vector defining which loci are transversions
#' @rdname loci_transversions
#' @export
#' @examples
#' example_gt <- example_gen_tibble()
#' example_gt %>% loci_transversions()
loci_transversions <- function(.x, ...) {
  UseMethod("loci_transversions", .x)
}

#' @export
#' @rdname loci_transversions
loci_transversions.tbl_df <- function(.x, ...) {
  # TODO this is a hack to deal with the class being dropped when going through
  # group_map
  stopifnot_gen_tibble(.x)
  check_allele_alphabet(.x$genotypes)
  loci_transversions(.x$genotypes, ...)
}


#' @export
#' @rdname loci_transversions
loci_transversions.vctrs_bigSNP <- function(.x, ...) {
  rlang::check_dots_empty()
  transversions <- function(loci_df) {
    (((loci_df$allele_ref == "A") & (loci_df$allele_alt == "T")) |
      ((loci_df$allele_ref == "T") & (loci_df$allele_alt == "A")) | # nolint start
      ((loci_df$allele_ref == "C") & (loci_df$allele_alt == "G")) |
      ((loci_df$allele_ref == "G") & (loci_df$allele_alt == "C"))) # nolint end
  }
  transversions(show_loci(.x))
}
