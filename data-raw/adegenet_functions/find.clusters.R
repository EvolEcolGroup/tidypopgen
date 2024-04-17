#' Run K-clustering on principal components
#'
#' This function implements the clustering procedure used in Discriminant
#' Analysis of Principal Components (DAPC, Jombart et al. 2010).
#' This procedure consists in running successive K-means with an
#' increasing number of clusters (k), after transforming data using
#' a principal component analysis (PCA). For each model,
#' a statistical measure of goodness of fit (by default, BIC)
#' is computed, which allows to choose the optimal k.
#' See details for a description of how to select the optimal k
#' and vignette("adegenet-dapc") for a tutorial.
#' @param x an object of class `gt_pca`, generated with [gt_pca()].
#' @param n_pca number of principal components to be fed to the LDA.
#' @param n_clust number of clusters
#' @param max_n_clust maximum number of clusters
#' @param n_iter number of iterations
#' @param n_start
#' @returns a [`gt_pca`] object with an additional element 'cluster',
#' which is a list with elements:
#' - 'method' the clustering method (either kmeans or ward)
#' - 'n_pca' number of principal components used for clustering
#' - 'k' the k values explored by the function
#' - 'WSS' within sum of squares for each k
#' - 'AIC' the AIC for each k
#' - 'BIC' the BIC for each k
#' - 'groups' a list, with each element giving the group assignments for a given k
#' @export

gt_pca_find_clusters <- function(x = NULL, n_pca = NULL,
                                 n_clusters = c(1,round(nrow(x$scores)/10)),
                                 method=c("kmeans","ward"),
                                   n_iter = 1e5, n_start = 10){

  if (is.null(x) | !inherits(x, "gt_pca")){
    stop("a 'gt_pca' object is required")
  }

  ## CHECKS ##
  method <- match.arg(method)

  ## SOME GENERAL VARIABLES ##
  N <- nrow(x$scores)
  # if we don't give the number of pca, then use them all
  if(is.null(n_pca)){
    n_pca <- length(anole_pca$eig)
    message("'n.pca' was not set: all ", n_pca, " components were used")
  }

  # unpack number of clusters
  if (length(n_clusters)==1){
    nbClust <- n_clusters
  } else if (length(n_clusters)==2){
    nbClust <- n_clusters[1]:n_clusters[2]
  } else {
    stop("'n_clusters' should be either a single value, or the minimum and maximum to be tested")
  }

  # set up list to store all results
  cluster_list <- list(method = method,
                       n_pca = n_pca,
                       k = nbClust,
                       WSS = numeric(0),
                       AIC = numeric(0),
                       BIC = numeric(0),
                       groups = list())

  # cluster of 1 is a special case
  if(nbClust[1]==1){
    nbClust = nbClust[-1]
    add_clust_of_1 <-TRUE
  } else {
    add_clust_of_1 <-FALSE
  }
  # vector to store within sum of squares
  WSS <- numeric(0)

  # TODO this is an obvious place where to parallelise
  for(i in 1:length(nbClust)){
    if (method == "kmeans") {
      ## kmeans clustering (original method)
      groups_assignments <- kmeans(x$scores[,seq_len(n_pca)],
                                   centers = nbClust[i],
                                   iter.max = n_iter, nstart = n_start)$cluster
    } else {
      ## ward clustering
      groups_assignments <- cutree(hclust(dist(x$scores[,seq_len(n_pca)])^2,
                                          method = "ward.D2"), k = nbClust[i])
    }
    WSS[i] <- .compute.wss(x$scores[,seq_len(n_pca)], groups_assignments)
    cluster_list$groups[[i]]<- groups_assignments

  }
  # compute statistics of goodnees of fit
  if (add_clust_of_1){
    WSS.ori <- sum(apply(x$scores[,seq_len(n_pca)], 2, function(v) sum((v-mean(v))^2) ))
    WSS <- c(WSS.ori,WSS)
    # add the classification for 1 cluster (they are all 1!)
    cluster_list$groups <- append(cluster_list$groups,
           list(setNames(rep(1,N),names(cluster_list$groups[[2]]))),
           after =0)
    nbClust <- c(1,nbClust)
  }
  cluster_list$AIC <- N*log(WSS/N) + 2 * nbClust
  cluster_list$BIC <- N*log(WSS/N) + log(N) *nbClust
  cluster_list$WSS <- WSS

  names(cluster_list$groups)<-nbClust

  x$clusters <- cluster_list
  class(x) <- c("gt_pca_clust", class(x))
  return(x)

}

## Compute within sum of squares from a matrix 'x' and a factor 'f'
.compute.wss <- function(x, f) {
  x.group.mean <- apply(x, 2, tapply, f, mean)
  sum((x - x.group.mean[as.character(f),])^2)
}



##TODO this function should add to the previous object
gt_pca_clust_choose_n <- function(x, stat = c("BIC", "AIC", "WSS"),
                                 criterion = c("diffNgroup", "min", "goesup",
                                                  "smoothNgoesup", "goodfit"),
                                 quiet=FALSE){
  if (!inherits(x, "gt_pca_clust")){
    stop("'x' should be a 'gt_pca_clusters' object generated with 'gt_pca_find_clusters()'")
  }

  stat <- match.arg(stat)
  criterion <- match.arg(criterion)

  if(criterion=="min") {
    n.clust <- which.min(x$clusters[[stat]])
  }
  if(criterion=="goesup") {
    ## temp <- diff(myStat)
    ## n.clust <- which.max( which( (temp-min(temp))<max(temp)/1e4))
    n.clust <- min(which(diff(x$clusters[[stat]])>0))
  }
  if(criterion=="goodfit") {
    temp <- min(x$clusters[[stat]]) + 0.1*(max(x$clusters[[stat]]) - min(x$clusters[[stat]]))
    n.clust <- min( which(x$clusters[[stat]] < temp))-1
  }
  if(criterion=="diffNgroup") {
    temp <- cutree(hclust(dist(diff(x$clusters[[stat]])), method="ward.D"), k=2)
    goodgrp <- which.min(tapply(diff(x$clusters[[stat]]), temp, mean))
    n.clust <- max(which(temp==goodgrp))+1
  }
  if(criterion=="smoothNgoesup") {
    temp <- x$clusters[[stat]]
    temp[2:(length(x$clusters[[stat]])-1)] <- sapply(1:(length(x$clusters[[stat]])-2),
                                         function(i) mean(x$clusters[[stat]][c(i,i+1,i+2)]))
    n.clust <- min(which(diff(temp)>0))
  }
  if(!quiet){
    message("Using ",stat, " with criterion ", criterion,": ", n.clust, " clusters")
  }
  x$n_clust <- n.clust
  return(x)
}
