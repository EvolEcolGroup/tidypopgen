#' Count the number of alleles or individuals per locus
#'
#' This function counts the number of alleles or individual for each locus
#' (excluding NAs for any given genotype).
#' It is not
#' meant to be used as a general summary, but it is used within other functions to
#' compute intermediate quantities.
#' When ploidy varies across individuals, the outputs of this function depend
#' on whether the information units are individuals, or
#' alleles within individuals (see details for [snpbin_list_sums()]).
#'
#' This function is a modified version
#' of [adegenet::glNA], recoded to work on lists of `SNPbin` objects as used
#' in the `genotypes` column of [gen_tibble].
#' @author Thibaut Jombart for the original [adegenet::glNA], modified
#' by Andrea Manica for 'tidypopgen'
#' @param .x a list of [`adegenet::SNPbin`] objects (usually the `genotype` column of
#' a [`gen_tibble`] object),
#' @param alleles_as_units a logical indicating whether alleles are considered
#' as units (i.e., a diploid genotype equals two samples, a triploid, three,
#' etc.) or whether individuals are considered as units of information.
#' @returns a vector of counts of alleles or individuals
#' @export


snpbin_list_n <- function(.x, alleles_as_units=TRUE){
  if (alleles_as_units) {
    return(sum(show_ploidy(.x)) -
             snpbin_list_count_na(.x,alleles_as_units = TRUE))

  } else {
    return(length(.x)-
      snpbin_list_count_na(.x,alleles_as_units = FALSE))
  }
}
