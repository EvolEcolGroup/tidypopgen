#' Estimate allele frequencies at each each locus
#'
#' Allele frequencies can be estimates as minimum allele frequencies (MAF) with
#' `loci_maf()` or the frequency of the alternate allele (with `loci_alt_freq()`).
#' The latter are in line with the genotypes matrix (e.g. as extracted by
#'  [`show_loci()`]). Most users will be in interested in the MAF, but the
#'  raw frequencies might be useful when computing aggregated statistics.
#'
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotypes` column of
#' a [`gen_tibble`] object),
#' or a [`gen_tibble`].
#' @param ... other arguments passed to specific methods, currently unused.
#' @returns a vector of frequencies, one per locus
#' @rdname loci_alt_freq
#' @export
loci_alt_freq <- function(.x, ...) {
  UseMethod("loci_alt_freq", .x)
}

#' @export
#' @rdname loci_alt_freq
loci_alt_freq.tbl_df <- function(.x, ...) {
  #TODO this is a hack to deal with the class being dropped when going through group_map
  stopifnot_gen_tibble(.x)
  loci_alt_freq(.x$genotypes, ...)
}


#' @export
#' @rdname loci_alt_freq
loci_alt_freq.vctrs_bigSNP <- function(.x, ...) {
  rlang::check_dots_empty()
  stopifnot_diploid(.x)
  # if we have diploid
  if (attr(.x,"ploidy")==2){
    loci_alt_freq_diploid(.x)
  }
}

#' @export
#' @rdname loci_alt_freq
loci_alt_freq.grouped_df <- function(.x, ...) {
  # TODO this is seriously inefficient, we need to cast it into a big_apply problem
  # of maybe it isn't that bad...
  group_map(.x, .f=~loci_alt_freq(.x,, ...))
}

#' @rdname loci_alt_freq
#' @export
loci_maf <- function(.x, ...) {
  UseMethod("loci_maf", .x)
}

#' @export
#' @rdname loci_alt_freq
loci_maf.tbl_df <- function(.x, ...) {
  #TODO this is a hack to deal with the class being dropped when going through group_map
  stopifnot_gen_tibble(.x)
  loci_maf(.x$genotypes, ...)
}

#' @export
#' @rdname loci_alt_freq
loci_maf.vctrs_bigSNP <- function(.x, ...) {
  freq <- loci_alt_freq(.x,)
  freq[freq>0.5 & !is.na(freq)] <- 1 - freq[freq>0.5 & !is.na(freq)]
  freq
}

#' @export
#' @rdname loci_alt_freq
loci_maf.grouped_df <- function(.x, ...) {
  # TODO this is seriously inefficient, we need to cast it into a big_apply problem
  # of maybe it isn't that bad...
  group_map(.x, .f=~loci_maf(.x, ...))
}

# function to estimate frequencies for diploid
loci_alt_freq_diploid <- function(.x){
  # get the FBM
  geno_fbm <- attr(.x,"bigsnp")$genotypes
  # rows (individuals) that we want to use
  rows_to_keep <- vctrs::vec_data(.x)
  # as long as we have more than one individual
  if (length(rows_to_keep)>1){
    col_counts <- bigstatsr::big_counts(geno_fbm,
                                        ind.row = rows_to_keep,
                                        ind.col = attr(.x,"loci")$big_index)
    means_from_counts <- function(x){
      (x[2]+x[3]*2)/((x[1]+x[2]+x[3])*2)
    }
    freq <- apply(col_counts, 2, means_from_counts)
  } else { # if we have a single individual
    freq <-geno_fbm[rows_to_keep,attr(.x,"loci")$big_index] /2
  }
  # sort out frequencies for diploids
  # @TODO get ploidy from the tibble once we have it
  # @TODO for mixed ploidy, we will need a different approach
  freq
}

loci_alt_freq_polyploid <- function(.x, ...){
  stop("not implemented yet")


}
