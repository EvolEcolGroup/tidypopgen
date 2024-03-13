#' Show the genotypes of a `gen_tibble`
#'
#' Extract the genotypes (as a matrix) from a  `gen_tibble`.
#' @param .x a list of [`adegenet::SNPbin`] objects (usually the `genotype` column of
#' a [`gen_tibble`] object),
#' or a [`gen_tibble`].
#' @param indiv_indices indices of individuals
#' @param loci_indices indices of loci
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
show_genotypes.tbl_df <- function(.x, indiv_indices=NULL, loci_indices=NULL, ...){
  stopifnot_gen_tibble(.x)
  # extract the column and hand it over to its method
  show_genotypes(.x$genotypes, indiv_indices=indiv_indices, loci_indices=loci_indices, ...)
}

#' @export
#' @rdname show_genotypes
show_genotypes.vctrs_bigSNP <- function(.x, indiv_indices=NULL, loci_indices=NULL, ...){
  rlang::check_dots_empty()
  if (is.null(indiv_indices) & is.null(loci_indices)){
    attr(.x,"bigsnp")$genotypes[,attr(.x,"loci")$big_index]
  } else if (is.null(indiv_indices)){
    attr(.x,"bigsnp")$genotypes[,attr(.x,"loci")$big_index[loci_indices]]
  } else if (is.null(loci_indices)){
    attr(.x,"bigsnp")$genotypes[indiv_indices,]
  } else {
    attr(.x,"bigsnp")$genotypes[indiv_indices,attr(.x,"loci")$big_index[loci_indices]]
  }
}
