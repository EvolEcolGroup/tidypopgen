#' Estimate allele frequencies at each each locus
#'
#' Estimate the frequency of the alternate allele at each locus.
#'
#' @param .x a list of [`adegenet::SNPbin`] objects (usually the `genotype` column of
#' a [`gen_tibble`] object),
#' or a [`gen_tibble`].
#' @param ... other arguments passed to specific methods.
#' @returns a vector of frequencies, one per locus
#' @rdname loci_freq
#' @export
loci_freq <- function(.x, ...) {
  UseMethod("loci_freq", .x)
}

#' @param alleles_as_units a logical indicating whether alleles are considered
#' as units (i.e., a diploid genotype equals two samples, a triploid, three,
#' etc.) or whether individuals are considered as units of information.
#' @param use_c a logical indicating whether compiled C code should be used
#' (TRUE) or not (FALSE, default).
#' @export
#' @rdname loci_freq
loci_freq.tbl_df <- function(.x, ..., minor = TRUE, alleles_as_units = TRUE, use_c = FALSE) {
  #TODO this is a hack to deal with the class being dropped when going through group_map
  stopifnot_gen_tibble(.x)
  loci_freq(.x$genotypes, ..., minor = minor, alleles_as_units = alleles_as_units, use_c = use_c)
}


#' @export
#' @rdname loci_freq
loci_freq.list <- function(.x, ..., minor = TRUE, alleles_as_units = TRUE, use_c = FALSE) {
  rlang::check_dots_empty()
  stopifnot_snpbin_list(.x)
  freq <- .genotypes_means(.x, alleles_as_units = alleles_as_units, use_c = use_c)
  if (minor){
    # TODO this fails because freq can be NA!!!
    freq[freq>0.5 & !is.na(freq)] <- 1 - freq[freq>0.5 & !is.na(freq)]
  }
  freq
}

#' @export
#' @rdname loci_freq
loci_freq.grouped_df <- function(.x, ..., minor = TRUE, alleles_as_units = TRUE, use_c = FALSE) {
  group_map(.x, .f=~loci_freq(.x, minor = minor, alleles_as_units = alleles_as_units, use_c = use_c))
}

