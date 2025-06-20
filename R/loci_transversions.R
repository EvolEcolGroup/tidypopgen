#' Find transversions
#'
#' Use the loci table to define which loci are transversions
#'
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotype` column of
#' a [`gen_tibble`] object),
#' or a [`gen_tibble`].
#' @param .col the column to be used when a tibble (or grouped tibble is passed
#' directly to the function). This defaults to "genotypes" and can only take
#' that value. There is no need for the user to set it, but it is included to
#' resolve certain tidyselect operations.
#' @param ... other arguments passed to specific methods.
#' @returns a logical vector defining which loci are transversions
#' @rdname loci_transversions
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#' example_gt %>% loci_transversions()
loci_transversions <- function(.x, .col = "genotypes", ...) {
  UseMethod("loci_transversions", .x)
}

#' @export
#' @rdname loci_transversions
loci_transversions.tbl_df <- function(.x, .col = "genotypes", ...) {
  stopifnot_gen_tibble(.x)
  .col <- rlang::enquo(.col) %>%
    rlang::quo_get_expr() %>%
    rlang::as_string()
  # confirm that .col is "genotypes"
  if (.col != "genotypes") {
    stop("loci_transversions only works with the genotypes column")
  }
  check_allele_alphabet(.x$genotypes)
  loci_transversions(.x$genotypes, ...)
}


#' @export
#' @rdname loci_transversions
loci_transversions.vctrs_bigSNP <- function(.x, .col = "genotypes", ...) {
  rlang::check_dots_empty()
  transversions <- function(loci_df) {
    (((loci_df$allele_ref == "A") & (loci_df$allele_alt == "T")) |
      ((loci_df$allele_ref == "T") & (loci_df$allele_alt == "A")) | # nolint start
      ((loci_df$allele_ref == "C") & (loci_df$allele_alt == "G")) |
      ((loci_df$allele_ref == "G") & (loci_df$allele_alt == "C"))) # nolint end
  }
  transversions(show_loci(.x))
}
