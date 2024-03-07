##########################
#' @method find.clusters gen_tibble
#' @export
##########################
find.clusters.gentibble <- function(x, clust = NULL, n.pca = NULL, n.clust = NULL,
                                   method = c("kmeans", "ward"),
                                   stat = c("BIC", "AIC", "WSS"),
                                   choose.n.clust = TRUE,
                                   criterion = c("diffNgroup", "min","goesup", "smoothNgoesup", "goodfit"),
                                   max.n.clust = round(nInd(x)/10), n.iter = 1e5, n.start = 10,
                                   scale = FALSE, pca.select = c("nbEig","percVar"),
                                   perc.pca = NULL, glPca = NULL, ...){

  ## CHECKS ##
  stat <- match.arg(stat)
  pca.select <- match.arg(pca.select)


  ## SOME GENERAL VARIABLES ##
  N <- nrow(x)
  min.n.clust <- 2


  ## PERFORM PCA ##
  REDUCEDIM <- is.null(glPca)

  if(REDUCEDIM){ # if no glPca provided
    maxRank <- min(c(nrow(x), length(show_loci_names(x))))
    pcaX <- gt_pca(x, center = TRUE, scale = scale, nf=maxRank, loadings=FALSE, returnDotProd = FALSE, ...)
  } else {
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


  ## convert PCA
  pcaX <- .glPca2dudi(pcaX)


  ## CALL DATA.FRAME METHOD
  res <- find.clusters(pcaX$li, clust=clust, n.pca=n.pca, n.clust=n.clust,
                       stat=stat, max.n.clust=max.n.clust, n.iter=n.iter, n.start=n.start,
                       choose.n.clust=choose.n.clust, method = method,
                       criterion=criterion, center=FALSE, scale=FALSE, dudi=pcaX)
  return(res)
} # end find.clusters.genlight



