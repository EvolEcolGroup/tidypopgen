##########################
#' @method find.clusters gen_tibble
#' @export
##########################

gt_find_clusters <- function(x, gt_pca=NULL, clust = NULL, n.pca = NULL, n.clust = NULL,
                                   method = c("kmeans", "ward"),
                                   max.n.clust = round(nrow(x)/10), n.iter = 1e5, n.start = 10,
                                   scale = FALSE){

  if (is.null(gt_pca) | !inherits(gt_pca, "gt_pca")){
    stop("a 'gt_pca' object is required")
  }

  ## CHECKS ##
  method <- match.arg(method)

  ## SOME GENERAL VARIABLES ##
  N <- nrow(gt_pca$scores)
  min.n.clust <- 2
  max.n.clust <- max(max.n.clust, 2)

  if(is.null(n.pca)){
    n.pca <- tidy(gt_pca ,matrix="eigenvalues") %>% dplyr::filter(.data$cumulative<0.9) %>% nrow() +1
    message("'n.pca' was not set: ", n.pca, " components were used")
  }

  nbClust <- min.n.clust:max.n.clust
  cluster_list <- list(k = seq_len(max.n.clust),
                       WSS = numeric(0),
                       AIC = numeric(0),
                       BIC = numeric(0),
                       groups = list())
  WSS <- numeric(0)

  for(i in 1:length(nbClust)){
    if (method == "kmeans") {
      ## kmeans clustering (original method)
      groups_assignments <- kmeans(gt_pca$scores[,seq_len(n.pca)], centers = nbClust[i], iter.max = n.iter, nstart = n.start)$cluster
    } else {
      ## ward clustering
      groups_assignments <- cutree(hclust(dist(gt_pca$scores[,seq_len(n.pca)])^2, method = "ward.D2"), k = nbClust[i])
    }
    WSS[i] <- .compute.wss(gt_pca$scores[,seq_len(n.pca)], groups_assignments)
    cluster_list$groups[[i+1]]<- groups_assignments

  }
  # compute statistics of goodnees of fit
  WSS.ori <- sum(apply(gt_pca$scores[,seq_len(n.pca)], 2, function(v) sum((v-mean(v))^2) ))
  cluster_list$AIC <- N*log(c(WSS.ori,WSS)/N) + 2*c(1,nbClust)
  cluster_list$BIC <- N*log(c(WSS.ori,WSS)/N) + log(N) *c(1,nbClust)
  cluster_list$WSS <- c(WSS.ori, WSS)
  names(cluster_list$groups)<-seq_len(1:max.n.clust)

  class(cluster_list) <- c("gt_clusters", class(cluster_list))
  return(cluster_list)

}




## Compute within sum of squares from a matrix 'x' and a factor 'f'
.compute.wss <- function(x, f) {
  x.group.mean <- apply(x, 2, tapply, f, mean)
  sum((x - x.group.mean[as.character(f),])^2)
}
