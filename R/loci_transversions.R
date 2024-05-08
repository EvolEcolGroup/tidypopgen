#' Estimates the sum of genotypes at each each locus
#'
#' Estimate the sum of the alternate allele at each locus. This is unlikely to be useful
#' directly, but it is used by other functions that compute various statistics.
#'
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotype` column of
#' a [`gen_tibble`] object),
#' or a [`gen_tibble`].
#' @param ... other arguments passed to specific methods.
#' @returns a vector of frequencies, one per locus
#' @rdname loci_transversions
#' @export
loci_transversions <- function(.x, ...) {
  UseMethod("loci_transversions", .x)
}

#' @export
#' @rdname loci_transversions
loci_transversions.tbl_df <- function(.x, ...) {
  #TODO this is a hack to deal with the class being dropped when going through group_map
  stopifnot_gen_tibble(.x)
  loci_transversions(.x$genotypes, ...)
}


#' @export
#' @rdname loci_transversions
loci_transversions.vctrs_bigSNP <- function(.x, ...) {
  rlang::check_dots_empty()
  transversions <- function(loci_df){
    (((loci_df$allele_ref=="A") & (loci_df$allele_alt=="T")) |
      ((loci_df$allele_ref=="T") & (loci_df$allele_alt=="A")) |
       ((loci_df$allele_ref=="C") & (loci_df$allele_alt=="G")) |
      ((loci_df$allele_ref=="G") & (loci_df$allele_alt=="C")))

  }
  transversions(show_loci(.x))
}

#' @export
#' @rdname loci_transversions
loci_transversions.grouped_df <- function(.x, ...) {
  group_map(.x, .f=~loci_transversions(.x))
}

