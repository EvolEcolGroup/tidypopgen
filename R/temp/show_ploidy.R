#' Show the loci information of a `gen_tibble`
#'
#' Extract the information on ploidy from a `gen_tibble`.
#' @param .x a list of [`adegenet::SNPbin`] objects (usually the `genotype` column of
#' a [`gen_tibble`] object),
#' or a [`gen_tibble`].
#' @param ... currently unused.
#' @returns a [`tibble::tibble`] of information (see[`gen_tibble`] for details
#' on compulsory columns that will always be present)
#' @rdname show_ploidy
#' @export
show_ploidy <- function(.x, ...) {
  UseMethod("show_ploidy", .x)
}

#' @export
#' @rdname show_ploidy
show_ploidy.gen_tbl <- function(.x, ...){
  stopifnot_gen_tibble(.x)
  show_ploidy(.x$genotypes, ...)
}

#' @export
#' @rdname show_ploidy
show_ploidy.list <- function(.x, ...){
  rlang::check_dots_empty()
  stopifnot_snpbin_list(.x)
  unlist(lapply(.x,adegenet::ploidy))
}
