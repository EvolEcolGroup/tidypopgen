#' Show the ploidy information of a `gen_tibble`
#'
#' Extract the ploidy information from a  `gen_tibble`.
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotype` column of
#' a [`gen_tibble`] object),
#' or a [`gen_tibble`].
#' @param ... currently unused.
#' @returns the ploidy (0 indicates mixed ploidy)
#' @rdname show_ploidy
#' @export
show_ploidy <- function(.x, ...) {
  UseMethod("show_ploidy", .x)
}

#' @export
#' @rdname show_ploidy
show_ploidy.tbl_df <- function(.x, ...){
  stopifnot_gen_tibble(.x)
  # extract the column and hand it over to its method
  show_ploidy(.x$genotypes)
}


#' @export
#' @rdname show_ploidy
show_ploidy.vctrs_bigSNP <- function(.x, ...){
  rlang::check_dots_empty()
  attr(.x,"ploidy")
}
