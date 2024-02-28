#' Show the loci information of a `gen_tibble`
#'
#' Extract the information on loci from a  `gen_tibble`.
#' @param .x a list of [`adegenet::SNPbin`] objects (usually the `genotype` column of
#' a [`gen_tibble`] object),
#' or a [`gen_tibble`].
#' @param ... currently unused.
#' @returns a [`tibble::tibble`] of information (see[`gen_tibble`] for details
#' on compulsory columns that will always be present)
#' @rdname show_loci
#' @export
show_loci <- function(.x, ...) {
  UseMethod("show_loci", .x)
}

#' @param .col if `.x` is a [`gen_tibble`], the column containing the genotypes
#' (usually `genotypes`)
#' @export
#' @rdname show_loci
show_loci.gen_tbl <- function(.x, .col, ...){
  # extract the column and hand it over to its method
  show_loci.list(.x[[rlang::ensym(.col)]], ...)
}

#' @export
#' @rdname show_loci
show_loci.list <- function(.x, ...){
  attr(.x,"loci")
}
