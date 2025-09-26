#' Show the loci information of a `gen_tibble`
#'
#' Extract and set the information on loci from a  `gen_tibble`.
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotype` column of
#' a [`gen_tibble`] object),
#' or a [`gen_tibble`].
#' @param ... currently unused.
#' @returns a [`tibble::tibble`] of information (see [`gen_tibble`] for details
#' on compulsory columns that will always be present)
#' @rdname show_loci
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' example_gt %>% show_loci()
show_loci <- function(.x, ...) {
  UseMethod("show_loci", .x)
}

#' @export
#' @rdname show_loci
show_loci.tbl_df <- function(.x, ...) {
  stopifnot_gen_tibble(.x)
  # extract the column and hand it over to its method
  show_loci(.x$genotypes, ...)
}


#' @export
#' @rdname show_loci
show_loci.vctrs_bigSNP <- function(.x, ...) {
  rlang::check_dots_empty()
  attr(.x, "loci")
}


#' @export
#' @param value a data.frame or tibble of loci information to replace the
#'   current one.
#' @rdname show_loci
"show_loci<-" <- function(.x, value) {
  UseMethod("show_loci<-", .x)
}

#' @export
#' @rdname show_loci
"show_loci<-.tbl_df" <- function(.x, value) {
  stopifnot_gen_tibble(.x)
  # extract the column and hand it over to its method
  show_loci(.x$genotypes) <- value
  .x
}

#' @export
#' @rdname show_loci
"show_loci<-.vctrs_bigSNP" <- function(.x, value) {
  # test for validity of loci
  value <- validate_loci(value)
  if (nrow(show_loci(.x)) != nrow(value)) {
    stop(paste(
      "the replacement loci tibble does not have the same number",
      "of loci as the original one"
    ))
  }
  attr(.x, "loci") <- tibble::as_tibble(value)
  .x
}
