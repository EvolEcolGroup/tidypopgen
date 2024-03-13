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
show_genotypes.tbl_df <- function(.x, ind_indices=NULL, loci_indices=NULL, ...){
  stopifnot_gen_tibble(.x)
  # extract the column and hand it over to its method
  show_genotypes(.x$genotypes, ind_indices=ind_indices, loci_indices=loci_indices, ...)
}

#' @export
#' @rdname show_genotypes
show_genotypes.vctrs_bigSNP <- function(.x, ind_indices=NULL, loci_indices=NULL, ...){
  rlang::check_dots_empty()
  if (is.null(ind_indices) & is.null(loci_indices)){
    attr(.x,"bigsnp")$genotypes[,attr(.x,"loci")$big_index]
  } else if (is.null(ind_indices)){
    attr(.x,"bigsnp")$genotypes[,attr(.x,"loci")$big_index[loci_indices]]
  } else if (is.null(loci_indices)){
    attr(.x,"bigsnp")$genotypes[ind_indices,]
  } else {
    attr(.x,"bigsnp")$genotypes[ind_indices,attr(.x,"loci")$big_index[loci_indices]]
  }
}
