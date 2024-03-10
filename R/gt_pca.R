#' Principal Component Analysis for `gen_tibble` object
#'
#' This function implements Principal Component Analysis for `gen_tibble`. It
#' is a modified version of [adegenet::glPca()], and returns an object that is
#' of class `glPca`, meaning that it is possible to use plotting and summary
#' functions from the
#' the package `adegenet` (for details on those function, see the help page
#' for [adegenet::glPca()].
#' @param .x	a [gen_tibble] object
#' @param center a logical indicating whether the numbers of alleles should
#' be centered; defaults to TRUE
#' @param scale	a logical indicating whether the numbers of alleles should
#' be scaled; defaults to FALSE
#' @param nf an integer indicating the number of principal components to be
#' retained; if NULL, all components with an eigen value larger than 1e-10 are kept.
#' @param loadings a logical indicating whether loadings of the alleles should
#' be computed (TRUE, default), or not (FALSE). Vectors of loadings are not
#' always useful, and can take a large amount of RAM when millions of SNPs are considered.
#' @param alleles_as_units a logical indicating whether alleles are considered as
#' units (i.e., a diploid genotype equals two samples, a triploid, three, etc.)
#' or whether individuals are considered as units of information.
#' @param use_c a logical indicating whether compiled C code should be used for
#' faster computations; this option cannot be used alongside parallel option.
#' @param parallel a logical indicating whether multiple cores -if available-
#' should be used for the computations (TRUE), or not (FALSE, default); requires
#' the package parallel to be installed (see details); this option cannot be
#' used alongside use_c option.
#' @param n_cores	if parallel is TRUE, the number of cores to be used in the
#' computations; if NULL, then the maximum number of cores available on the computer is used.
#' @param return_dot_prod	a logical indicating whether the matrix of dot products
#' between individuals should be returned (TRUE) or not (FALSE, default).
#' @param mat_dot_prod an optional matrix of dot products between individuals,
#' NULL by default. This option is used internally to speed up computation
#' time when re-running the same PCA several times. Leave this argument as NULL
#' unless you really know what you are doing.
#' @returns an object of class `gt_pca` (a subclass of [`adegenet::glPca`], with the following components:
#' call	the matched call.
#' eig a numeric vector of eigenvalues.
#' scores a matrix of principal components, containing the coordinates of each individual (in row) on each principal axis (in column).
#' loadings (optional) a matrix of loadings, containing the loadings of each SNP (in row) for each principal axis (in column).
#' @export

gt_pca <- function(.x, center=TRUE, scale=FALSE, nf=NULL, loadings=TRUE, alleles_as_units=FALSE,
                  use_c=TRUE, parallel=FALSE, n_cores=NULL,
                  return_dot_prod=FALSE, mat_dot_prod=NULL){
  #browser()
  nInd <- nrow(.x)
  nLoc <- nrow(show_loci(.x))
  ## COMPUTE MEANS AND VARIANCES ##
  if(center) {
    vecMeans <- snpbin_list_means(.x$genotypes, alleles_as_units=alleles_as_units)
    if(any(is.na(vecMeans))) stop("NAs detected in the vector of means")
  }

  if(scale){
    vecVar <- snpbin_list_vars(.x$genotypes, alleles_as_units=alleles_as_units)
    if(any(is.na(vecVar))) stop("NAs detected in the vector of variances")
  }


  myPloidy <- show_ploidy(.x)

  ## NEED TO COMPUTE DOT PRODUCTS ##
  if(is.null(mat_dot_prod)){

    ## == if non-C code is used ==
    if(!use_c){
      ## if(parallel && !require(parallel)) stop("parallel package requested but not installed")
      if(parallel && is.null(n_cores)){
        n_cores <- parallel::detectCores()
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
      allComb <- utils::combn(1:nInd, 2)
      if(parallel){
        allProd <- unlist(parallel::mclapply(1:ncol(allComb), function(i) dotProd(.x[[allComb[1,i]]], .x$genotypes[[allComb[2,i]]], myPloidy[allComb[1,i]], myPloidy[allComb[2,i]]),
                                             mc.cores=n_cores, mc.silent=TRUE, mc.cleanup=TRUE, mc.preschedule=FALSE))
      } else {
        allProd <- unlist(lapply(1:ncol(allComb), function(i) dotProd(.x$genotypes[[allComb[1,i]]], .x$genotypes[[allComb[2,i]]], myPloidy[allComb[1,i]], myPloidy[allComb[2,i]]) ))
      }
      allProd <- allProd / nInd # assume uniform weights

      ## shape result as a matrix
      attr(allProd,"Size") <- nInd
      attr(allProd,"Diag") <- FALSE
      attr(allProd,"Upper") <- FALSE
      class(allProd) <- "dist"
      allProd <- as.matrix(allProd)

      ## compute the diagonal
      if(parallel){
        temp <- unlist(parallel::mclapply(1:nInd, function(i) dotProd(.x$genotypes[[i]], .x$genotypes[[i]], myPloidy[i], myPloidy[i]),
                                          mc.cores=n_cores, mc.silent=TRUE, mc.cleanup=TRUE, mc.preschedule=FALSE))/nInd
      } else {
        temp <- unlist(lapply(1:nInd, function(i) dotProd(.x$genotypes[[i]], .x$genotypes[[i]], myPloidy[i], myPloidy[i]) ))/nInd
      }
      diag(allProd) <- temp
    } else { # === use C computations ====
      allProd <- snpbin_list_dot_prod(.x$genotypes, center=center, scale=scale, alleles_as_units=alleles_as_units, parallel=parallel, n_cores=n_cores)/nInd
    }
  } else { # END NEED TO COMPUTE DOTPROD
    if(!all(dim(mat_dot_prod)==nInd)) stop("mat_dot_prod has wrong dimensions.")
    allProd <- mat_dot_prod
  }

  ## PERFORM THE ANALYSIS ##
  ## eigenanalysis
  eigRes <- eigen(allProd, symmetric=TRUE, only.values=FALSE)
  rank <- sum(eigRes$values > 1e-12)
  eigRes$values <- eigRes$values[1:rank]
  eigRes$vectors <- eigRes$vectors[, 1:rank, drop=FALSE]

  if (is.null(nf)){
    nf <- rank
  }

  ## rescale PCs
  res <- list()
  res$eig <- eigRes$values
  nf <- min(nf, sum(res$eig>1e-10))
  ##res$matprod <- allProd # for debugging

  ## use: li = XQU = V\Lambda^(1/2)
  eigRes$vectors <- eigRes$vectors * sqrt(nInd) # D-normalize vectors
  res$scores <- sweep(eigRes$vectors[, 1:nf, drop=FALSE],2, sqrt(eigRes$values[1:nf]), FUN="*")


  ## GET LOADINGS ##
  ## need to decompose X^TDV into a sum of n matrices of dim p*r
  ## but only two such matrices are represented at a time
  if(loadings){
    if(scale) {
      vecSd <- sqrt(vecVar)
    }
    res$loadings <- matrix(0, nrow=nLoc, ncol=nf) # create empty matrix
    ## use: c1 = X^TDV
    ## and X^TV = A_1 + ... + A_n
    ## with A_k = X_[k-]^T v[k-]
    for(k in 1:nInd){
      temp <- as.integer(.x$genotypes[[k]]) / myPloidy[k]
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

    res$loadings <- res$loadings / nInd # don't forget the /n of X_tDV
    res$loadings <- sweep(res$loadings, 2, sqrt(eigRes$values[1:nf]), FUN="/")
  }

  #browser()
  ## FORMAT OUTPUT ##
  colnames(res$scores) <- paste("PC", 1:nf, sep="")
  if(!is.null(.x$id)){
    rownames(res$scores) <- .x$id
  } else {
    rownames(res$scores) <- 1:nInd
  }

  if(!is.null(res$loadings)){
    colnames(res$loadings) <- paste("Axis", 1:nf, sep="")
    # TODO quick hack, bring back the labels AM
    # if(!is.null(show_loci_names(x_gt)) & !is.null(alleles(x))){
    #   rownames(res$loadings) <- paste(show_loci_names(x_gt),alleles(x), sep=".")
    # } else {
      rownames(res$loadings) <- 1:nLoc
    # }
  }

  if(return_dot_prod){
    res$dotProd <- allProd
    rownames(res$dotProd) <- colnames(res$dotProd) <- .x$id
  }

  # save center and scale information
  if(center) {
    res$cent <- vecMeans
  }

  if(scale){
    res$norm <- sqrt(vecVar)
  }

  res$call <- match.call()

  class(res) <- c("gt_pca","glPca")

  return(res)

}

# a print method
#' @method print gt_pca
#' @export
print.gt_pca <- function(x, ...){
  cat(" === PCA of gen_tibble object ===")
  cat("\nClass: list of type gt_pca glPca")
  cat("\nCall ($call):")
  print(x$call)
  cat("\nEigenvalues ($eig):\n", round(utils::head(x$eig,6),3), ifelse(length(x$eig)>6, "...\n", "\n") )
  cat("\nPrincipal components ($scores):\n matrix with", nrow(x$scores), "rows (individuals) and", ncol(x$scores), "columns (axes)", "\n")
  if(!is.null(x$loadings)){
    cat("\nPrincipal axes ($loadings):\n matrix with", nrow(x$loadings), "rows (SNPs) and", ncol(x$loadings), "columns (axes)", "\n")
  }
  if(!is.null(x$dotProd)){
    cat("\nDot products between individuals ($dotProd):\n matrix with", nrow(x$dotProd), "rows and", ncol(x$dotProd), "columns", "\n")
  }
  cat("\n")
}
