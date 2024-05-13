#' Find transitions
#'
#' Use the loci table to define which loci are transitions
#'
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotype` column of
#' a [`gen_tibble`] object),
#' or a [`gen_tibble`].
#' @param ... other arguments passed to specific methods.
#' @returns a vector of frequencies, one per locus
#' @rdname loci_transitions
#' @export
loci_transitions <- function(.x, ...) {
  UseMethod("loci_transitions", .x)
}

#' @export
#' @rdname loci_transitions
loci_transitions.tbl_df <- function(.x, ...) {
  #TODO this is a hack to deal with the class being dropped when going through group_map
  stopifnot_gen_tibble(.x)
  loci_transitions(.x$genotypes, ...)
}


#' @export
#' @rdname loci_transitions
loci_transitions.vctrs_bigSNP <- function(.x, ...) {
  rlang::check_dots_empty()

  !loci_transversions(.x)
}

#' @export
#' @rdname loci_transitions
loci_transitions.grouped_df <- function(.x, ...) {
  group_map(.x, .f=~loci_transitions(.x))
}

