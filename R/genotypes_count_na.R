#' Count the number of missing genotypes
#'
#' This function counts the number of missing genotypes for each SNP.
#' When ploidy varies across individuals, the outputs of this function depend
#' on whether the information units are individuals, or
#' alleles within individuals (see details).
#'
#' @param .x a list of [`adegenet::SNPbin`] objects (usually the `genotype` column of
#' a [`gen_tibble`] object),
#' @param alleles_as_unit a logical indicating whether alleles are considered
#' as units (i.e., a diploid genotype equals two samples, a triploid, three,
#' etc.) or whether individuals are considered as units of information.
#' @returns a vector of counts of NAs
#' @export


genotypes_count_na <- function(.x, alleles_as_unit=TRUE){

  naPosi <- lapply(.x,adegenet::NA.posi)
  nLoci <- nrow(attr(.x,"loci"))
  ploidy <- attr(.x,"ploidy")

  ## DEFAULT, VECTOR-WISE PROCEDURE ##
  res <- integer(nLoci)
  temp <- naPosi

  ## NAs in allele sampling
  if(alleles_as_unit){
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
