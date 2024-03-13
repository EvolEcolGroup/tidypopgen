#' Compute the means of genotypes
#'
#' This function computes the means of genotypes.
#'
#' This function is a modified version
#' of [adegenet::glMean], recoded to work on lists of `SNPbin` objects as used
#' in the `genotypes` column of [gen_tibble].
#' @author Thibaut Jombart for the original [adegenet::glSum], modified
#' by Andrea Manica for 'tidypopgen'
#' @param .x a list of [`adegenet::SNPbin`] objects (usually the `genotype` column of
#' a [`gen_tibble`] object),
#' @param alleles_as_units a logical indicating whether alleles are considered
#' as units (i.e., a diploid genotype equals two samples, a triploid, three,
#' etc.) or whether individuals are considered as units of information.
#' @param use_c a logical indicating whether compiled C code should be used
#' (TRUE) or not (FALSE, default).
#' @returns a vector of means
#' @export


snpbin_list_means <- function(.x, alleles_as_units = TRUE, use_c = FALSE){

  ## DEFAULT, VECTOR-WISE PROCEDURE ##
  if(alleles_as_units){ # use alleles
    ploidy <- show_ploidy(.x)
    N <- sum(ploidy) - snpbin_list_count_na(.x, alleles_as_units=TRUE)
  } else { # use relative frequencies of individuals
    nInd <- length(.x)
    N <- nInd - snpbin_list_count_na(.x, alleles_as_units=FALSE)
  }
  res <- snpbin_list_sums(.x, alleles_as_units=alleles_as_units, use_c = use_c)/N

  names(res) <- attr(.x,"loci")$name
  return(res)
}
