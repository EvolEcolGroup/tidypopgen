##########################
#' @method find.clusters gen_tibble
#' @export
##########################

gt_find_clusters <- function(x, gt_pca=NULL, clust = NULL, n.pca = NULL, n.clust = NULL,
                                   method = c("kmeans", "ward"),
                                   stat = c("BIC", "AIC", "WSS"),
                                   choose.n.clust = TRUE,
                                   criterion = c("diffNgroup", "min","goesup", "smoothNgoesup", "goodfit"),
                                   max.n.clust = round(nrow(x)/10), n.iter = 1e5, n.start = 10,
                                   scale = FALSE){

  if (is.null(gt_pca) | !inherits(gt_pca, "gt_pca")){
    stop("a 'gt_pca' object is required")
  }

  ## CHECKS ##
  stat <- match.arg(stat)


  ## SOME GENERAL VARIABLES ##
  N <- nrow(x)
  min.n.clust <- 2

  if(is.null(n.pca)){
    n.pca <- tidy(gt_pca ,matrix="eigenvalues") %>% dplyr::filter(.data$cumulative<0.9) %>% nrow() +1
  }

  ## convert PCA
  pcaX <- adegenet:::.glPca2dudi(gt_pca)


  ## CALL DATA.FRAME METHOD
  res <- find.clusters(pcaX$li, clust=clust, n.pca=n.pca, n.clust=n.clust,
                       stat=stat, max.n.clust=max.n.clust, n.iter=n.iter, n.start=n.start,
                       choose.n.clust=choose.n.clust, method = method,
                       criterion=criterion, center=FALSE, scale=FALSE, dudi=pcaX)
  return(res)
} # end find.clusters.genlight



