#################
## dapc.genlight
#################
#' @method dapc gen_tibble
#' @export
dapc.gen_tibble <- function(x, pop=NULL, n.pca=NULL, n.da=NULL,
                          scale=FALSE,  var.contrib=TRUE, var.loadings=FALSE, pca.info=TRUE,
                          pca.select=c("nbEig","percVar"), perc.pca=NULL, glPca=NULL, ...){
  pca.select <- match.arg(pca.select)

  if(is.null(pop)) {
    pop.fac <- x$population
  } else {
    pop.fac <- pop
  }

  if(is.null(pop.fac)) stop("x does not include pre-defined populations, and `pop' is not provided")



  ## PERFORM PCA ##
  REDUCEDIM <- is.null(glPca)

  if(REDUCEDIM){ # if no glPca provided
    maxRank <- min(c(nrow(x), length(show_loci_names(x))))
    pcaX <- gt_pca(x, center = TRUE, scale = scale, nf=maxRank, loadings=FALSE, returnDotProd = TRUE, ...)
  }

  if(!REDUCEDIM){ # else use the provided glPca object
    if(is.null(glPca$loadings) & var.contrib) {
      warning("Contribution of variables requested but glPca object provided without loadings.")
      var.contrib <- FALSE
    }
    pcaX <- glPca
  }

  if(is.null(n.pca)){
    cumVar <- 100 * cumsum(pcaX$eig)/sum(pcaX$eig)
  }


  ## select the number of retained PC for PCA
  if(!REDUCEDIM){
    myCol <- rep(c("black", "lightgrey"), c(ncol(pcaX$scores),length(pcaX$eig)))
  } else {
    myCol <- "black"
  }

  if(is.null(n.pca) & pca.select=="nbEig"){
    plot(cumVar, xlab="Number of retained PCs", ylab="Cumulative variance (%)", main="Variance explained by PCA", col=myCol)
    cat("Choose the number PCs to retain (>=1): ")
    n.pca <- as.integer(readLines(con = getOption('adegenet.testcon'), n = 1))
  }

  if(is.null(perc.pca) & pca.select=="percVar"){
    plot(cumVar, xlab="Number of retained PCs", ylab="Cumulative variance (%)", main="Variance explained by PCA", col=myCol)
    cat("Choose the percentage of variance to retain (0-100): ")
    nperc.pca <- as.numeric(readLines(con = getOption('adegenet.testcon'), n = 1))
  }

  ## get n.pca from the % of variance to conserve
  if(!is.null(perc.pca)){
    n.pca <- min(which(cumVar >= perc.pca))
    if(perc.pca > 99.999) n.pca <- length(pcaX$eig)
    if(n.pca<1) n.pca <- 1
  }

  if(!REDUCEDIM){
    if(n.pca > ncol(pcaX$scores)) {
      n.pca <- ncol(pcaX$scores)
    }
  }


  ## recompute PCA with loadings if needed
  if(REDUCEDIM){
    pcaX <- gt_pca(x, center = TRUE, scale = scale, nf=n.pca, loadings=var.contrib, matDotProd = pcaX$dotProd)
  }


  ## keep relevant PCs - stored in XU
  N <- nrow(x)
  X.rank <- sum(pcaX$eig > 1e-14)
  n.pca <- min(X.rank, n.pca)
  if(n.pca >= N) n.pca <- N-1

  U <- pcaX$loadings[, 1:n.pca, drop=FALSE] # principal axes
  XU <- pcaX$scores[, 1:n.pca, drop=FALSE] # principal components
  XU.lambda <- sum(pcaX$eig[1:n.pca])/sum(pcaX$eig) # sum of retained eigenvalues
  names(U) <- paste("PCA-pa", 1:ncol(U), sep=".")
  names(XU) <- paste("PCA-pc", 1:ncol(XU), sep=".")


  ## PERFORM DA ##
  ldaX <- lda(XU, pop.fac, tol=1e-30) # tol=1e-30 is a kludge, but a safe (?) one to avoid fancy rescaling by lda.default
  lda.dim <- sum(ldaX$svd^2 > 1e-10)
  ldaX$svd <- ldaX$svd[1:lda.dim]
  ldaX$scaling <- ldaX$scaling[,1:lda.dim,drop=FALSE]

  if(is.null(n.da)){
    barplot(ldaX$svd^2, xlab="Linear Discriminants", ylab="F-statistic", main="Discriminant analysis eigenvalues", col=heat.colors(length(levels(pop.fac))) )
    cat("Choose the number discriminant functions to retain (>=1): ")
    n.da <- as.integer(readLines(con = getOption('adegenet.testcon'), n = 1))
  }

  n.da <- min(n.da, length(levels(pop.fac))-1, n.pca, sum(ldaX$svd>1e-10)) # can't be more than K-1 disc. func., or more than n.pca
  n.da <- round(n.da)
  predX <- predict(ldaX, dimen=n.da)


  ## BUILD RESULT
  res <- list()
  res$n.pca <- n.pca
  res$n.da <- n.da
  res$tab <- XU
  res$grp <- pop.fac
  res$var <- XU.lambda
  res$eig <- ldaX$svd^2
  res$loadings <- ldaX$scaling[, 1:n.da, drop=FALSE]
  res$means <- ldaX$means
  res$ind.coord <-predX$x
  res$grp.coord <- apply(res$ind.coord, 2, tapply, pop.fac, mean)
  res$prior <- ldaX$prior
  res$posterior <- predX$posterior
  res$assign <- predX$class
  res$call <- match.call()


  ## optional: store loadings of variables
  if(pca.info){
    res$pca.loadings <- as.matrix(U)
    res$pca.cent <- snpbin_list_means(x,alleles_as_units =FALSE)
    if(scale) {
      res$pca.norm <- sqrt(snpbin_list_vars(x,alleles_as_units = FALSE))
    } else {
      res$pca.norm <- rep(1, nLoc(x))
    }
    res$pca.eig <- pcaX$eig
  }

  ## optional: get loadings of variables
  if(var.contrib || var.loadings){
    var.load <- as.matrix(U) %*% as.matrix(ldaX$scaling[,1:n.da,drop=FALSE])

    if(var.contrib){
      f1 <- function(x){
        temp <- sum(x*x)
        if(temp < 1e-12) return(rep(0, length(x)))
        return(x*x / temp)
      }
      res$var.contr <- apply(var.load, 2, f1)
    }
    if(var.loadings) res$var.load <- var.load
  }

  class(res) <- "dapc"
  return(res)
} # end dapc.genlight

