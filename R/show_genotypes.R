#' Show the genotypes of a `gen_tibble`
#'
#' Extract the genotypes (as a matrix) from a  `gen_tibble`.
#' @param x a [`gen_tibble`] object.
#' @param ... currently unused
#' @returns a [`tibble::tibble`] of information (see[`gen_tibble`] for details
#' on compulsory columns that will always be present)
#' @rdname show_genotypes
#' @export
show_genotypes <- function(x, ...) {
  UseMethod("show_genotypes", x)
}

#' @export
#' @rdname show_genotypes
show_genotypes.gen_tbl <- function(x, ...){
  res <- unlist(lapply(x$genotypes, as.integer))
  res <- matrix(res, ncol=nrow(attr(x,"loci")), nrow = nrow(x), byrow=TRUE)
  colnames(res) <- attr(x,"loci")$name
  rownames(res) <- x$id
  res
}

#' @export
#' @rdname show_genotypes
show_genotypes.list <- function(x, ...){
  if (!inherits(x[[1]],"SNPBin")){ # for the sake of speed, we only check the first element
    stop("x is not a list of SNPBin objects")
  }
  res <- unlist(lapply(x, as.integer))
  res <- matrix(res, nrow = length(x), byrow=TRUE)
  res
}
