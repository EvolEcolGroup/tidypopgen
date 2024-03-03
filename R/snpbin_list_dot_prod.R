#' Compute the the dot products of loci
#'
#' This function compute the dot products between centere/scaled vectors of
#' loci. It is not
#' meant to be used as a general summary, but it is used within other functions to
#' compute intermediate quantities.
#' When ploidy varies across individuals, the outputs of this function depend
#' on whether the information units are individuals, or
#' alleles within individuals (see details for [snpbin_list_sums()]).
#'
#' This function is a modified version
#' of [adegenet::glDotProd], recoded to work on lists of `SNPbin` objects as used
#' in the `genotypes` column of [gen_tibble].
#' @author Thibaut Jombart for the original [adegenet::glDotProd], modified
#' by Andrea Manica for 'tidypopgen'
#' @param .x a list of [`adegenet::SNPbin`] objects (usually the `genotype` column of
#' a [`gen_tibble`] object).
#' @param center a logical indicating whether SNPs should be centred to mean zero.
#' @param scale a logical indicating whether SNPs should be scaled to unit variance.
#' @param alleles_as_units a logical indicating whether alleles are considered
#' as units (i.e., a diploid genotype equals two samples, a triploid, three,
#' etc.) or whether individuals are considered as units of information.
#' @param parallel a logical indicating whether multiple cores -if
#' available- should be used for the computations (TRUE, default), or
#' not (FALSE); requires the package `parallel` to be installed
#' (see details); this option cannot be used alongside use_c option
#' @param n_cores if `parallel` is TRUE, the number of cores to
#' be used in the computations; if NULL, then the maximum number of
#' cores available on the computer is used.

#' @returns a vector of counts of NAs
#' @export


snpbin_list_dot_prod <- function(.x, center=FALSE, scale=FALSE, alleles_as_units=FALSE,
                        parallel=FALSE, n_cores=NULL){

    if(parallel && is.null(n_cores)){
      n_cores <- parallel::detectCores()
    }


    ## STORE USEFUL INFO ##
    N <- length(.x)
    nInd <- N
    nLoci <- nrow(attr(.x,"loci"))
#    ind.names <- indNames(x)


    if(!parallel){ # DO NOT USE MULTIPLE CORES
      ## GET INPUTS TO C PROCEDURE ##
      if(center){
        mu <- snpbin_list_means(.x,alleles_as_units=alleles_as_units)
      } else {
        mu <- rep(0, nLoci)
      }

      if(scale){
        s <- sqrt(snpbin_list_vars(.x,alleles_as_units=alleles_as_units))
        if(any(s<1e-10)) {
          warning("Null variances have been detected; corresponding alleles won't be standardized.")
        }
      } else {
        s <- rep(1, nLoci)
      }

      vecbyte <- unlist(lapply(.x, function(e) e$snp))
      nbVec <- sapply(.x, function(e) length(e$snp))
      naPosi <- lapply(.x,adegenet::NA.posi)
      nbNa <- sapply(naPosi, length)
      naPosi <- unlist(naPosi)
      lowerTriSize <- (nInd*(nInd-1))/2
      resSize <- lowerTriSize + nInd

#      return(list(vecbyte, nbVec, length(.x[[1]]@snp[[1]]), nbNa, naPosi, nInd, nLoci, show_ploidy(.x),
#                  as.double(mu), as.double(s), as.integer(!alleles_as_units), double(resSize)))
      ## CALL C FUNCTION ##return
      temp <- .C("GLdotProd", vecbyte, nbVec, length(.x[[1]]@snp[[1]]), nbNa, naPosi, nInd, nLoci, show_ploidy(.x),
                 as.double(mu), as.double(s), as.integer(!alleles_as_units), double(resSize), PACKAGE="adegenet")[[12]]
#      return(temp)
    } else { # USE MULTIPLE CORES
      stop("this has not been implemented yet!")
      # @FIXME we need to implement something equivalent to seploc!!!!
      # I don't think this is efficient, as it creates a bit copy in memory
      # we would be better off indexing the blocks
      .x <- seploc(.x, n.block = n_cores) # one block per core (x is now a list of genlight)
      temp <- list()
      i <- 0
      for(block in .x){
        i <- i+1
        ## GET INPUTS TO C PROCEDURE ##
        if(center){
          mu <- snpbin_list_means(block,alleles_as_units=alleles_as_units)
        } else {
          mu <- rep(0, nrow(attr(block,"loci")))
        }

        if(scale){
          s <- sqrt(snpbin_list_vars(block,alleles_as_units=alleles_as_units))
          if(any(s<1e-10)) {
            warning("Null variances have been detected; corresponding alleles won't be standardized.")
          }
        } else {
          s <- rep(1, nrow(attr(block,"loci")))
        }

        vecbyte <- unlist(lapply(block, function(e) e$snp))
        nbVec <- sapply(block, function(e) length(e$snp))
        naPosi <- lapply(block,adegenet::NA.posi)
        nbNa <- sapply(naPosi, length)
        naPosi <- unlist(naPosi)

        nInd_block <- length(block)
        lowerTriSize <- (nInd_block*(nInd_block-1))/2
        resSize <- lowerTriSize + nInd_block

        ## CALL C FUNCTION ##
        temp[[i]] <- .C("GLdotProd", vecbyte, nbVec, length(block[[1]]@snp[[1]]), nbNa, naPosi, nInd_block,
                        nrow(attr(block,"loci")), attr(block,"ploidy"),
                        as.double(mu), as.double(s), as.integer(!alleles_as_units), double(resSize), PACKAGE="adegenet")[[12]]
      }


      ## POOL BLOCK RESULTS TOGETHER ##
      temp <- Reduce("+", temp)
    }

    res <- temp[1:lowerTriSize]
    attr(res,"Size") <- N
    attr(res,"Diag") <- FALSE
    attr(res,"Upper") <- FALSE
    class(res) <- "dist"
    res <- as.matrix(res)
    diag(res) <- temp[(lowerTriSize+1):length(temp)]

    return(res)
  }


