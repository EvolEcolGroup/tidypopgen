#' Estimate individual heterozygosity
#'
#' Estimate heterozygosity for each individual.
#'
#' @param .x a list of [`adegenet::SNPbin`] objects (usually the `genotype` column of
#' a [`gen_tibble`] object),
#' or a [`gen_tibble`].
#' @param ... currently unused.
#' @returns a vector of heterozygosities, one per individuals in the [`gen_tibble`]
#' @rdname heterozygosity
#' @export
heterozygosity <- function(.x, ...) {
  UseMethod("heterozygosity", .x)
}

#' @param .col if `.x` is a [`gen_tibble`], the column containing the genotypes
#' (usually `genotypes`)
#' @export
#' @rdname heterozygosity
heterozygosity.gen_tbl <- function(.x, .col, ...){
  # extract the column and hand it over to its method
  heterozygosity(.x[[rlang::ensym(.col)]], ...)
}

#' @export
#' @rdname heterozygosity
heterozygosity.list <- function(.x, ...){
  if (!inherits(.x[[1]],"SNPbin")){ # for the sake of speed, we only check the first element
    stop(".x is not a list of SNPbin objects")
  }
  ## TODO we should implement loci means directly on the SNPbin objects
  rowMeans(show_genotypes(.x) == 1, na.rm = TRUE)
}
