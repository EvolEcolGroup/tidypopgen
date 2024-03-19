#' Test Hardy-Weinberg equilibrium at each locus
#'
#' Return the p-value from an exact test of HWE. It uses direclty the C++ code
#' from PLINK as provided in [HardyWeinberg::HWExactStats()]).
#'
#' NOTE There are no tests for this function yet! Unit tests are needed.
#'
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotypes` column of
#' a [`gen_tibble`] object),
#' or a [`gen_tibble`].
#' @param ... further arguments to pass to [HardyWeinberg::HWExact()].
#' @returns a vector of probabilities from HWE exact test, one per locus
#' @author based on code originally written by Jan Graffleman, optimised by Andrea Manica for speed
#' @rdname loci_hwe
#' @export
loci_hwe <- function(.x, ...) {
  UseMethod("loci_hwe", .x)
}

# We should write a cpp counts function. We can't use the snp_* family of functions
# as they ignore NAs

#' @export
#' @rdname loci_hwe
loci_hwe.tbl_df <- function(.x, ...) {
  #TODO this is a hack to deal with the class being dropped when going through group_map
  stopifnot_gen_tibble(.x)
  loci_hwe(.x$genotypes, ...)
}


#' @export
#' @rdname loci_hwe
loci_hwe.vctrs_bigSNP <- function(.x, ...) {
  #rlang::check_dots_empty()
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
      apply(geno_counts,2,HWExact_geno_col)
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
  group_map(.x, .f=~loci_hwe(.x))
}

HWExact_geno_col <- function(x){
  # use C++ unexported function from HardyWeinberg
  # it would be even better to do the counting in C++
  HardyWeinberg:::SNPHWE2(x[2],
                          x[1],
                          x[3],
                          midp=TRUE)
}
