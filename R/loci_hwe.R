#' Test Hardy-Weinberg equilibrium at each locus
#'
#' Return the p-value from an exact test of HWE.
#'
#' This function uses the original C++ algorithm
#' from PLINK 1.90.
#'
#' NOTE There are no tests for this function yet! Unit tests are needed.
#'
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotypes` column of
#' a [`gen_tibble`] object),
#' or a [`gen_tibble`].
#' @param mid_p boolean on whether the mid-p value should be computed.
#' Default is TRUE, as in PLINK.
#' @param ... not used.
#' @returns a vector of probabilities from HWE exact test, one per locus
#' @author the C++ algorithm was written by Christopher Chang for PLINK 1.90, based on
#' original code by Jan Wigginton (the code was released under GPL3).
#' @rdname loci_hwe
#' @export
loci_hwe <- function(.x, ...) {
  UseMethod("loci_hwe", .x)
}


#' @export
#' @rdname loci_hwe
loci_hwe.tbl_df <- function(.x, mid_p = TRUE, ...) {
  #TODO this is a hack to deal with the class being dropped when going through group_map
  stopifnot_gen_tibble(.x)
  loci_hwe(.x$genotypes, mid_p = mid_p, ...)
}


#' @export
#' @rdname loci_hwe
loci_hwe.vctrs_bigSNP <- function(.x, mid_p = TRUE, ...) {
  rlang::check_dots_empty()
  stopifnot_diploid(.x)
  # get the FBM
  geno_fbm <- attr(.x,"bigsnp")$genotypes
  # rows (individuals) that we want to use
  rows_to_keep <- vctrs::vec_data(.x)
  # as long as we have more than one individual
  if (length(rows_to_keep)>1){
    # col means for submatrix (all rows, only some columns)
    colHWE_sub <- function(X, ind, rows_to_keep) {
      #apply(X[rows_to_keep, ind], 2, HWExact_geno_vec)
      geno_counts <- bigstatsr::big_counts(X,ind.row = rows_to_keep, ind.col = ind)
      apply(geno_counts,2,HWExact_geno_col, mid_p = mid_p)
    }
    hwe_p <- bigstatsr::big_apply(geno_fbm, a.FUN = colHWE_sub,
                                 rows_to_keep = rows_to_keep,
                                 ind=attr(.x,"loci")$big_index,
                                 a.combine = 'c')
  } else { # if we have a single individual
    stop ("Not implemented for a single individual")
  }
    hwe_p
}

#' @export
#' @rdname loci_hwe
loci_hwe.grouped_df <- function(.x, ...) {
  # TODO this is seriously inefficient, we need to cast it into a big_apply problem
  # of maybe it isn't that bad...
  group_map(.x, .f=~loci_hwe(.x, mid_p = mid_p, ...))
}

HWExact_geno_col <- function(x, mid_p){
  # it would be even better to use it direclty in a C function that does the counting
  SNPHWE2(x[2],
                          x[1],
                          x[3],
                          midp=mid_p)
}
