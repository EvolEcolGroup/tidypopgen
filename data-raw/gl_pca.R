########
## glPca
########
##
## PCA for genlight objects
##
glPca <- function(x, center=TRUE, scale=FALSE, nf=NULL, loadings=TRUE, alleles_as_units=FALSE,
                  useC=TRUE, parallel=FALSE, n.cores=NULL,
                  returnDotProd=FALSE, matDotProd=NULL){

  ## COMPUTE MEANS AND VARIANCES ##
  if(center) {
    vecMeans <- .genotype_means(x, alleles_as_units=alleles_as_units)
    if(any(is.na(vecMeans))) stop("NAs detected in the vector of means")
  }

  if(scale){
    vecVar <- .genotype_vars(x, alleles_as_units=alleles_as_units)
    if(any(is.na(vecVar))) stop("NAs detected in the vector of variances")
  }


  myPloidy <- show_ploidy(x)

  ## NEED TO COMPUTE DOT PRODUCTS ##
  if(is.null(matDotProd)){

    ## == if non-C code is used ==
    if(!useC){
      ## if(parallel && !require(parallel)) stop("parallel package requested but not installed")
      if(parallel && is.null(n.cores)){
        n.cores <- parallel::detectCores()
      }


      ## COMPUTE DOT PRODUCTS BETWEEN GENOTYPES ##
      ## to be fast, a particular function is defined for each case of centring/scaling

      ## NO CENTRING / NO SCALING
      if(!center & !scale){
        dotProd <- function(a,b, ploid.a, ploid.b){ # a and b are two SNPbin objects
          a <- as.integer(a) / ploid.a
          a[is.na(a)] <- 0
          b <- as.integer(b) / ploid.b
          b[is.na(b)] <- 0
          return(sum( a*b, na.rm=TRUE))
        }
      }

      ## CENTRING / NO SCALING
      if(center & !scale){
        dotProd <- function(a,b, ploid.a, ploid.b){ # a and b are two SNPbin objects
          a <- as.integer(a) / ploid.a
          a[is.na(a)] <- vecMeans[is.na(a)]
          b <- as.integer(b) / ploid.b
          b[is.na(b)] <- vecMeans[is.na(b)]
          return(sum( (a-vecMeans) * (b-vecMeans), na.rm=TRUE) )
        }
      }


      ## NO CENTRING / SCALING (odd option...)
      if(!center & scale){
        dotProd <- function(a,b, ploid.a, ploid.b){ # a and b are two SNPbin objects
          a <- as.integer(a) / ploid.a
          a[is.na(a)] <- 0
          b <- as.integer(b) / ploid.b
          b[is.na(b)] <- 0
          return(sum( (a*b)/vecVar, na.rm=TRUE))
        }
      }


      ## CENTRING / SCALING
      if(center & scale){
        dotProd <- function(a,b, ploid.a, ploid.b){ # a and b are two SNPbin objects
          a <- as.integer(a) / ploid.a
          a[is.na(a)] <- vecMeans[is.na(a)]
          b <- as.integer(b) / ploid.b
          b[is.na(b)] <- vecMeans[is.na(b)]
          return( sum( ((a-vecMeans)*(b-vecMeans))/vecVar, na.rm=TRUE ) )
        }
      }


      ## COMPUTE ALL POSSIBLE DOT PRODUCTS (XX^T / n) ##
      allComb <- combn(1:nInd(x), 2)
      if(parallel){
        allProd <- unlist(parallel::mclapply(1:ncol(allComb), function(i) dotProd(x[[allComb[1,i]]], x$genotypes[[allComb[2,i]]], myPloidy[allComb[1,i]], myPloidy[allComb[2,i]]),
                                             mc.cores=n.cores, mc.silent=TRUE, mc.cleanup=TRUE, mc.preschedule=FALSE))
      } else {
        allProd <- unlist(lapply(1:ncol(allComb), function(i) dotProd(x$genotypes[[allComb[1,i]]], x$genotypes[[allComb[2,i]]], myPloidy[allComb[1,i]], myPloidy[allComb[2,i]]) ))
      }
      allProd <- allProd / nInd(x) # assume uniform weights

      ## shape result as a matrix
      attr(allProd,"Size") <- nInd(x)
      attr(allProd,"Diag") <- FALSE
      attr(allProd,"Upper") <- FALSE
      class(allProd) <- "dist"
      allProd <- as.matrix(allProd)

      ## compute the diagonal
      if(parallel){
        temp <- unlist(parallel::mclapply(1:nInd(x), function(i) dotProd(x$genotypes[[i]], x$genotypes[[i]], myPloidy[i], myPloidy[i]),
                                          mc.cores=n.cores, mc.silent=TRUE, mc.cleanup=TRUE, mc.preschedule=FALSE))/nInd(x)
      } else {
        temp <- unlist(lapply(1:nInd(x), function(i) dotProd(x$genotypes[[i]], x$genotypes[[i]], myPloidy[i], myPloidy[i]) ))/nInd(x)
      }
      diag(allProd) <- temp
    } else { # === use C computations ====
      allProd <- glDotProd(x, center=center, scale=scale, alleles_as_units=alleles_as_units, parallel=parallel, n.cores=n.cores)/nInd(x)
    }
  } else { # END NEED TO COMPUTE DOTPROD
    if(!all(dim(matDotProd)==nInd(x))) stop("matDotProd has wrong dimensions.")
    allProd <- matDotProd
  }

  ## PERFORM THE ANALYSIS ##
  ## eigenanalysis
  eigRes <- eigen(allProd, symmetric=TRUE, only.values=FALSE)
  rank <- sum(eigRes$values > 1e-12)
  eigRes$values <- eigRes$values[1:rank]
  eigRes$vectors <- eigRes$vectors[, 1:rank, drop=FALSE]

  ## scan nb of axes retained
  if(is.null(nf)){
    barplot(eigRes$values, main="Eigenvalues", col=heat.colors(rank))
    cat("Select the number of axes: ")
    nf <- as.integer(readLines(con = getOption('adegenet.testcon'), n = 1))
  }

  ## rescale PCs
  res <- list()
  res$eig <- eigRes$values
  nf <- min(nf, sum(res$eig>1e-10))
  ##res$matprod <- allProd # for debugging

  ## use: li = XQU = V\Lambda^(1/2)
  eigRes$vectors <- eigRes$vectors * sqrt(nInd(x)) # D-normalize vectors
  res$scores <- sweep(eigRes$vectors[, 1:nf, drop=FALSE],2, sqrt(eigRes$values[1:nf]), FUN="*")


  ## GET LOADINGS ##
  ## need to decompose X^TDV into a sum of n matrices of dim p*r
  ## but only two such matrices are represented at a time
  if(loadings){
    if(scale) {
      vecSd <- sqrt(vecVar)
    }
    res$loadings <- matrix(0, nrow=nLoc(x), ncol=nf) # create empty matrix
    ## use: c1 = X^TDV
    ## and X^TV = A_1 + ... + A_n
    ## with A_k = X_[k-]^T v[k-]
    for(k in 1:nInd(x)){
      temp <- as.integer(x$genotypes[[k]]) / myPloidy[k]
      if(center) {
        temp[is.na(temp)] <- vecMeans[is.na(temp)]
        temp <- temp - vecMeans
      } else {
        temp[is.na(temp)] <- 0
      }
      if(scale){
        temp <- temp/vecSd
      }

      res$loadings <- res$loadings + matrix(temp) %*% eigRes$vectors[k, 1:nf, drop=FALSE]
    }

    res$loadings <- res$loadings / nInd(x) # don't forget the /n of X_tDV
    res$loadings <- sweep(res$loadings, 2, sqrt(eigRes$values[1:nf]), FUN="/")
  }

  browser()
  ## FORMAT OUTPUT ##
  colnames(res$scores) <- paste("PC", 1:nf, sep="")
  if(!is.null(indNames(x))){
    rownames(res$scores) <- indNames(x)
  } else {
    rownames(res$scores) <- 1:nInd(x)
  }

  if(!is.null(res$loadings)){
    colnames(res$loadings) <- paste("Axis", 1:nf, sep="")
    if(!is.null(locNames(x)) & !is.null(alleles(x))){
      rownames(res$loadings) <- paste(locNames(x),alleles(x), sep=".")
    } else {
      rownames(res$loadings) <- 1:nLoc(x)
    }
  }

  if(returnDotProd){
    res$dotProd <- allProd
    rownames(res$dotProd) <- colnames(res$dotProd) <- indNames(x)
  }

  res$call <- match.call()

  class(res) <- "glPca"

  return(res)

} # glPca
