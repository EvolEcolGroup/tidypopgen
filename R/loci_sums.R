#' Estimates the sum of genotypes at each each locus
#'
#' Estimate the sum of the alternate allele at each locus.
#'
#' @param .x a list of [`adegenet::SNPbin`] objects (usually the `genotype` column of
#' a [`gen_tibble`] object),
#' or a [`gen_tibble`].
#' @param ... other arguments passed to specific methods.
#' @returns a vector of frequencies, one per locus
#' @rdname loci_sums
#' @export
loci_sums <- function(.x, ...) {
  UseMethod("loci_sums", .x)
}

#' @param minor a logical indicating whether we should give the frequencies of
#' the minor allele (TRUE, the default). If FALSE, the frequencies of the
#' alternate allele are given.
#' @export
#' @rdname loci_sums
loci_sums.tbl_df <- function(.x, ...) {
  #TODO this is a hack to deal with the class being dropped when going through group_map
  stopifnot_gen_tibble(.x)
  loci_sums(.x$genotypes, ..., minor = minor)
}


#' @export
#' @rdname loci_sums
loci_sums.vctrs_bigSNP <- function(.x, ...) {
  rlang::check_dots_empty()
  # get the FBM
  geno_fbm <- attr(.x,"bigsnp")$genotypes
  # rows (individuals) that we want to use
  rows_to_keep <- vctrs::vec_data(.x)
  # as long as we have more than one individual
  if (length(rows_to_keep)>1){
    # col means for submatrix (all rows, only some columns)
    colMeans_sub <- function(X, ind, rows_to_keep) {
      colSums(X[rows_to_keep, ind], na.rm=TRUE)
    }
    freq <- bigstatsr::big_apply(geno_fbm, a.FUN = colMeans_sub,
                                 rows_to_keep = rows_to_keep,
                                 ind=attr(.x,"loci")$big_index,
                                 a.combine = 'c')
  } else { # if we have a single individual
    freq <-geno_fbm[rows_to_keep,attr(.x,"loci")$big_index]
  }
  freq
}

#' @export
#' @rdname loci_sums
loci_sums.grouped_df <- function(.x, ...) {
  # TODO this is seriously inefficient, we need to cast it into a big_apply problem
  # of maybe it isn't that bad...
  group_map(.x, .f=~loci_sums(.x, minor = minor))
}

