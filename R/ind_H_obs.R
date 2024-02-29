#' Estimate individual observed heterozygosity
#'
#' Estimate heterozygosity for each individual.
#'
#' @param .x a list of [`adegenet::SNPbin`] objects (usually the `genotype` column of
#' a [`gen_tibble`] object),
#' or a [`gen_tibble`].
#' @param ... currently unused.
#' @returns a vector of heterozygosities, one per individuals in the [`gen_tibble`]
#' @rdname ind_H_obs
#' @export
ind_H_obs <- function(.x, ...) {
  UseMethod("ind_H_obs", .x)
}

#' @param .col if `.x` is a [`gen_tibble`], the column containing the genotypes
#' (usually `genotypes`)
#' @export
#' @rdname ind_H_obs
ind_H_obs.gen_tbl <- function(.x, .col, ...){
  # extract the column and hand it over to its method
  ind_H_obs(.x[[rlang::ensym(.col)]], ...)
}

#' @export
#' @rdname ind_H_obs
ind_H_obs.list <- function(.x, ...){
  if (!inherits(.x[[1]],"SNPbin")){ # for the sake of speed, we only check the first element
    stop(".x is not a list of SNPbin objects")
  }
  ## function to compute observed homozygosity in a SNPbin
  h_obs <- function(this_snpbin){
    mean(as.integer(this_snpbin)==1, na.rm = TRUE)
  }
  unlist(lapply(.x,h_obs))
}
