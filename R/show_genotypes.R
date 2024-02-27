#' Show the genotypes of a `gen_tibble`
#'
#' Extract the genotypes (as a matrix) from a  `gen_tibble`.
#' @param .data a [`gen_tibble`] object.
#' @param ... currently unused
#' @returns a [`tibble::tibble`] of information (see[`gen_tibble`] for details
#' on compulsory columns that will always be present)
#' @rdname show_genotypes
#' @export
show_genotypes <- function(.data, ...) {
  UseMethod("show_genotypes", .data)
}

#' @export
#' @rdname show_genotypes
show_genotypes.gen_tbl <- function(.data, ...){
  res <- unlist(lapply(.data$genotypes, as.integer))
  res <- matrix(res, ncol=nrow(attr(.data,"loci")), nrow = nrow(.data), byrow=TRUE)
  colnames(res) <- attr(.data,"loci")$name
  rownames(res) <- .data$id
  res
}

#' @export
#' @rdname show_genotypes
show_genotypes.list <- function(.data, ...){
  if (!inherits(.data[[1]],"SNPbin")){ # for the sake of speed, we only check the first element
    stop("x is not a list of SNPbin objects")
  }
  res <- unlist(lapply(.data, as.integer))
  res <- matrix(res, nrow = length(.data), byrow=TRUE)
  res
}
