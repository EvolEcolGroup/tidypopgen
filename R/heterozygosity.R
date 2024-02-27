#' Estimate individual heterozygosity
#'
#' Estimate heterozygosity for each individual.
#'
#' @param .data a [`gen_tibble`] object.
#' @param ... currently unused
#' @returns a vector of heterozygosities
#' @rdname heterozygosity
#' @export
heterozygosity <- function(.data, ...) {
  UseMethod("heterozygosity", .data)
}

#' @export
#' @rdname heterozygosity
heterozygosity.gen_tbl <- function(.data, ...){
  ## TODO we should implement loci means directly on the SNPbin objects
  rowMeans(show_genotypes(.data) == 1, na.rm = TRUE)
}

#' @export
#' @rdname heterozygosity
heterozygosity.list <- function(.data, ...){
  if (!inherits(.data[[1]],"SNPbin")){ # for the sake of speed, we only check the first element
    stop(".data is not a list of SNPbin objects")
  }
  ## TODO we should implement loci means directly on the SNPbin objects
  rowMeans(show_genotypes(.data) == 1, na.rm = TRUE)
}
