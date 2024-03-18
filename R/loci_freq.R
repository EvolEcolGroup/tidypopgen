#' Estimate allele frequencies at each each locus
#'
#' Estimate the frequency of the alternate allele at each locus.
#'
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotypes` column of
#' a [`gen_tibble`] object),
#' or a [`gen_tibble`].
#' @param ... other arguments passed to specific methods.
#' @returns a vector of frequencies, one per locus
#' @rdname loci_freq
#' @export
loci_freq <- function(.x, ...) {
  UseMethod("loci_freq", .x)
}

#' @param minor a logical indicating whether we should give the frequencies of
#' the minor allele (TRUE, the default). If FALSE, the frequencies of the
#' alternate allele are given.
#' @export
#' @rdname loci_freq
loci_freq.tbl_df <- function(.x, ..., minor = TRUE) {
  #TODO this is a hack to deal with the class being dropped when going through group_map
  stopifnot_gen_tibble(.x)
  loci_freq(.x$genotypes, ..., minor = minor)
}


#' @export
#' @rdname loci_freq
loci_freq.vctrs_bigSNP <- function(.x, ..., minor = TRUE) {
  rlang::check_dots_empty()
  # get the FBM
  geno_fbm <- attr(.x,"bigsnp")$genotypes
  # rows (individuals) that we want to use
  rows_to_keep <- vctrs::vec_data(.x)
  # as long as we have more than one individual
  if (length(rows_to_keep)>1){
    # col means for submatrix (all rows, only some columns)
    colMeans_sub <- function(X, ind, rows_to_keep) {
      colMeans(X[rows_to_keep, ind], na.rm=TRUE)
    }
    freq <- bigstatsr::big_apply(geno_fbm, a.FUN = colMeans_sub,
                                 rows_to_keep = rows_to_keep,
                                 ind=attr(.x,"loci")$big_index,
                                 a.combine = 'c')
  } else { # if we have a single individual
    freq <-geno_fbm[rows_to_keep,attr(.x,"loci")$big_index]
  }
  # sort out frequencies for diploids
  freq <- freq/2
  if (minor){
    freq[freq>0.5 & !is.na(freq)] <- 1 - freq[freq>0.5 & !is.na(freq)]
  }
  freq
}

#' @export
#' @rdname loci_freq
loci_freq.grouped_df <- function(.x, ..., minor = TRUE) {
  # TODO this is seriously inefficient, we need to cast it into a big_apply problem
  # of maybe it isn't that bad...
  group_map(.x, .f=~loci_freq(.x, minor = minor))
}

