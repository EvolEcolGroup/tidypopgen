#' Sum of genotypes
#'
#' This function sums the number of alternate alleles in each SNP. It is not
#' meant to be used as a general summary, but it is used within other functions to
#' compute intermediate quantities.
#' When ploidy varies across individuals, the outputs of this function depend
#' on whether the information units are individuals, or
#' alleles within individuals (see details).
#'
#' === On the unit of information ===
#' In the cases where individuals can have different ploidy, computation of
#' sums, means, etc. of allelic data depends on what we consider as a
#' unit of information.
#' To estimate e.g. allele frequencies, unit of information can be
#' considered as the allele, so that a diploid genotype contains two
#' samples, a triploid individual, three samples, etc. In such a case,
#' all computations are done directly on the number of alleles. This
#' corresponds to allele_as_unit = TRUE.
#'
#' However, when the focus is put on studying differences/similarities
#' between individuals, the unit of information is the individual, and all
#' genotypes possess the same information no matter what their ploidy is.
#' In this case, computations are made after standardizing individual
#' genotypes to relative allele frequencies. This corresponds
#' to allele_as_units = FALSE.
#'
#' Note that when all individuals have the same ploidy, this distinction
#' does not hold any more.
#'
#' This function is a modified version
#' of [adegenet::glSum], recoded to work on lists of `SNPbin` objects as used
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
#' @returns a [`tibble::tibble`] of information (see[`gen_tibble`] for details
#' on compulsory columns that will always be present)
#' @export

.genotypes_sums <- function (.x, alleles_as_units=TRUE, use_c=FALSE){

  nInd <- length(.x)
  nLoci <- nrow(attr(.x,"loci"))
  ploidy <- attr(.x,"ploidy")

  if(use_c){
      vecbyte <- unlist(lapply(.x, function(e) e$snp))
      nbVec <- sapply(.x, function(e) length(e$snp))
      naPosi <- lapply(.x,adegenet::NA.posi)
      nbNa <- sapply(naPosi, length)
      naPosi <- unlist(naPosi)
      if(alleles_as_units){
        ## use ploidy (sum absolute frequencies)
        res <- .C("GLsumInt", vecbyte, nbVec, length(.x[[1]]@snp[[1]]), nbNa, naPosi,
                nInd, nLoci, ploidy,
                integer(nLoci), PACKAGE="adegenet")[[9]]
      } else {
        ## sum relative frequencies
        res <- .C("GLsumFreq", vecbyte, nbVec, length(.x[[1]]@snp[[1]]), nbNa, naPosi,
                  nInd, nLoci, ploidy,
                  double(nLoci), PACKAGE="adegenet")[[9]]
      }
  } else {
    ## use ploidy (sum absolute frequencies)
    if(alleles_as_units){
      res <- integer(nLoci)
      for(e in .x){
        temp <- as.integer(e)
        temp[is.na(temp)] <- 0L
        res <- res + temp
      }
    } else {
      ## sum relative frequencies
      res <- numeric(nLoci)
      for(i in 1:nInd){
        temp <- as.integer(.x[[i]]) / ploidy[i]
        temp[is.na(temp)] <- 0
        res <- res + temp
      }
    }

  }
  names(res) <- attr(.x,"loci")$name
  return(res)

}
