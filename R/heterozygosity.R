#' Estimate individual heterozygosity
#'
#' Estimate heterozygosity for each individual.
#'
#' @param x a [`gen_tibble`] object.
#' @param ... currently unused
#' @returns a vector of heterozygosities
#' @rdname heterozygosity
#' @export
heterozygosity <- function(x, ...) {
  UseMethod("heterozygosity", x)
}

#' @export
#' @rdname heterozygosity
heterozygosity.gen_tbl <- function(x, ...){
  ## TODO we should implement loci means directly on the SNPbin objects
  rowMeans(show_genotypes(x) == 1, na.rm = TRUE)
}

#' @export
#' @rdname heterozygosity
heterozygosity.list <- function(x, ...){
  if (!inherits(x[[1]],"SNPBin")){ # for the sake of speed, we only check the first element
    stop("x is not a list of SNPBin objects")
  }
  ## TODO we should implement loci means directly on the SNPbin objects
  rowMeans(show_genotypes(x) == 1, na.rm = TRUE)
}
