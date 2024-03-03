#' Count the number of missing genotypes
#'
#' This function counts the number of missing genotypes for each locus. It is not
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
#' @returns a vector of counts of NAs
#' @export


snpbin_list_count_na <- function(.x, alleles_as_units=TRUE){

  naPosi <- lapply(.x,adegenet::NA.posi)
  nLoci <- nrow(attr(.x,"loci"))
  ploidy <- show_ploidy(.x)

  ## DEFAULT, VECTOR-WISE PROCEDURE ##
  res <- integer(nLoci)
  temp <- naPosi

  ## NAs in allele sampling
  if(alleles_as_units){
    for(i in 1:length(temp)){
      if(length(temp[[i]])>0){
        res[temp[[i]]] <- res[temp[[i]]] + ploidy[i]
      }
    }
  } else { ## NAs amongst individuals
    for(e in temp){
      if(length(e)>0){
        res[e] <- res[e] + 1
      }
    }
  }

  names(res) <- attr(.x,"loci")$name
  return(res)

}
