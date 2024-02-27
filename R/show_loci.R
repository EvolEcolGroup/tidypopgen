#' Show the loci information of a `gen_tibble`
#'
#' Extract the information on loci from a  `gen_tibble`.
#' @param .data a [`gen_tibble`] object.
#' @param ... currently unused
#' @returns a [`tibble::tibble`] of information (see[`gen_tibble`] for details
#' on compulsory columns that will always be present)
#' @rdname show_loci
#' @export
show_loci <- function(.data, ...) {
  UseMethod("show_loci", .data)
}

#' @export
#' @rdname show_loci
show_loci.gen_tbl <- function(.data, ...){
  attr(.data,"loci")
}
