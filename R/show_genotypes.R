#' Show the genotypes of a `gen_tibble`
#'
#' Extract the genotypes (as a matrix) from a  `gen_tibble`.
#' @param .x a list of [`adegenet::SNPbin`] objects (usually the `genotype` column of
#' a [`gen_tibble`] object),
#' or a [`gen_tibble`].
#' @param ... currently unused.
#' @returns a matrix of counts of the alternative alleles (see [show_loci()]) to
#' extract information on the alleles for those loci from a [`gen_tibble`].
#' @rdname show_genotypes
#' @export
show_genotypes <- function(.x, ...) {
  UseMethod("show_genotypes", .x)
}

#' @export
#' @rdname show_genotypes
show_genotypes.gen_tbl <- function(.x, ...){
  # extract the column and hand it over to its method
  #show_genotypes.list(.x[[rlang::ensym(.col)]])
  show_genotypes(.x$genotypes)
}

#' @export
#' @rdname show_genotypes
show_genotypes.list <- function(.x, ...){
  rlang::check_dots_empty()
  if (!inherits(.x[[1]],"SNPbin")){ # for the sake of speed, we only check the first element
    stop("x is not a list of SNPbin objects")
  }
  res <- unlist(lapply(.x, as.integer))
  res <- matrix(res, ncol=nrow(attr(.x,"loci")), nrow = length(.x), byrow=TRUE)
  colnames(res) <- attr(.x,"loci")$name
  res
}
