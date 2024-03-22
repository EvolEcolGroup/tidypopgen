#' Compute the variances of loci
#'
#' This function compute the variance for each locus (taking into account
#' missing genotypes). It is not
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
#' @param .x a vector of class `vctrs_bigSNP` (usually the `genotype` column of
#' a [`gen_tibble`] object),
#' @param alleles_as_units a logical indicating whether alleles are considered
#' as units (i.e., a diploid genotype equals two samples, a triploid, three,
#' etc.) or whether individuals are considered as units of information.
#' @param use_c a logical indicating whether compiled C code should be used
#' (TRUE) or not (FALSE, default).
#' @returns a vector of counts of NAs
#' @export


snpbin_list_vars <-
  function(.x,
           alleles_as_units = TRUE,
           use_c = FALSE) {
    nInd <- length(.x)
    nLoci <- nrow(attr(.x, "loci"))
    ploidy <- show_ploidy(.x)

    ## DEFAULT, VECTOR-WISE PROCEDURE ##
    res <- numeric(nLoci)

    if (alleles_as_units) {
      # use alleles
      N <-
        sum(ploidy) - snpbin_list_count_na(.x, alleles_as_units = TRUE)
      xbar <-
        snpbin_list_means(.x, alleles_as_units = TRUE, use_c = use_c)
      for (i in 1:nInd) {
        temp <- (as.integer(.x[[i]]) / ploidy[i] - xbar) ^ 2
        temp[is.na(temp)] <- 0
        res <- res + temp * ploidy[i]
      }
      res <- res / N
    } else {
      # use relative frequencies of individuals
      N <- nInd - snpbin_list_count_na(.x, alleles_as_units = FALSE)
      xbar <-
        snpbin_list_means(.x, alleles_as_units = FALSE, use_c = use_c)

      for (i in 1:nInd) {
        temp <- (as.integer(.x[[i]]) / ploidy[i] - xbar) ^ 2
        temp[is.na(temp)] <- 0L
        res <- res + temp
      }
      res <- res / N
    }

    names(res) <- attr(.x, "loci")$name
    return(res)

  }
