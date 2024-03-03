glDotProd2 <- function(x, center=FALSE, scale=FALSE, alleleAsUnit=FALSE,
                      parallel=FALSE, n.cores=NULL){
  if(!inherits(x, "genlight")) stop("x is not a genlight object")

  ## SOME CHECKS ##
  ## if(parallel && !require(parallel)) stop("parallel package requested but not installed")
  if(parallel && is.null(n.cores)){
    n.cores <- parallel::detectCores()
  }


  ## STORE USEFUL INFO ##
  N <- nInd(x)
  ind.names <- indNames(x)


  if(!parallel){ # DO NOT USE MULTIPLE CORES
    ## GET INPUTS TO C PROCEDURE ##
    if(center){
      mu <- glMean(x,alleleAsUnit=alleleAsUnit)
    } else {
      mu <- rep(0, nLoc(x))
    }

    if(scale){
      s <- sqrt(glVar(x,alleleAsUnit=alleleAsUnit))
      if(any(s<1e-10)) {
        warning("Null variances have been detected; corresponding alleles won't be standardized.")
      }
    } else {
      s <- rep(1, nLoc(x))
    }

    vecbyte <- unlist(lapply(x@gen, function(e) e$snp))
    nbVec <- sapply(x@gen, function(e) length(e$snp))
    nbNa <- sapply(NA.posi(x), length)
    naPosi <- unlist(NA.posi(x))
    lowerTriSize <- (nInd(x)*(nInd(x)-1))/2
    resSize <- lowerTriSize + nInd(x)

   return(list(vecbyte, nbVec, length(x@gen[[1]]@snp[[1]]), nbNa, naPosi, nInd(x), nLoc(x), ploidy(x),
               as.double(mu), as.double(s), as.integer(!alleleAsUnit), double(resSize)))
    ## CALL C FUNCTION ##
    temp <- .C("GLdotProd", vecbyte, nbVec, length(x@gen[[1]]@snp[[1]]), nbNa, naPosi, nInd(x), nLoc(x), ploidy(x),
               as.double(mu), as.double(s), as.integer(!alleleAsUnit), double(resSize), PACKAGE="adegenet")[[12]]
#    return(temp)
  } else { # USE MULTIPLE CORES
    x <- seploc(x, n.block = n.cores) # one block per core (x is now a list of genlight)
    temp <- list()
    i <- 0
    for(block in x){
      i <- i+1
      ## GET INPUTS TO C PROCEDURE ##
      if(center){
        mu <- glMean(block,alleleAsUnit=alleleAsUnit)
      } else {
        mu <- rep(0, nLoc(block))
      }

      if(scale){
        s <- sqrt(glVar(block,alleleAsUnit=alleleAsUnit))
        if(any(s<1e-10)) {
          warning("Null variances have been detected; corresponding alleles won't be standardized.")
        }
      } else {
        s <- rep(1, nLoc(block))
      }

      vecbyte <- unlist(lapply(block@gen, function(e) e$snp))
      nbVec <- sapply(block@gen, function(e) length(e$snp))
      nbNa <- sapply(NA.posi(block), length)
      naPosi <- unlist(NA.posi(block))
      lowerTriSize <- (nInd(block)*(nInd(block)-1))/2
      resSize <- lowerTriSize + nInd(block)

      ## CALL C FUNCTION ##
      temp[[i]] <- .C("GLdotProd", vecbyte, nbVec, length(block@gen[[1]]@snp[[1]]), nbNa, naPosi, nInd(block), nLoc(block), ploidy(block),
                      as.double(mu), as.double(s), as.integer(!alleleAsUnit), double(resSize), PACKAGE="adegenet")[[12]]
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

  colnames(res) <- rownames(res) <- ind.names
  return(res)
} # end glDotProd

