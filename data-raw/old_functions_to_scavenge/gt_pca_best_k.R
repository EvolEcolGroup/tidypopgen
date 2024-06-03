#' Find the best number of clusters based on principal components
#'
#' This function selects the best k value based on a chosen metric and criterion.
#' It is equivalent to plotting the metric against the *k* values, and selecting
#' the *k* that fulfils a given criterion (see details for an explanation of
#' each criterion). This function simply adds an element 'best_k' to the
#' `gt_cluster_pca` returned by [gt_cluster_pca()]. The choice can be
#' over-ridden simply by assigning a different value to that element (e.g.
#' for an object x and a desired *k* of 8, simply use `x$best_k <- 8`)
#'
#' The analysis of data simulated under various population genetics models
#' (see reference) suggested an ad-hoc rule for the selection of the optimal
#' number of clusters. First important result is that BIC seems more
#' efficient than AIC and WSS to select the appropriate number of clusters
#' (see example). The rule of thumb consists in increasing K until it
#' no longer leads to an appreciable improvement of fit (i.e., to a decrease
#' of BIC). In the most simple models (island models), BIC decreases until
#' it reaches the optimal K, and then increases. In these cases, the best rule
#' amounts to choosing the lowest K. In other models such as stepping
#' stones, the decrease of BIC often continues after the optimal K,
#' but is much less steep, so a change in slope can be taken as an indication
#' of where the best *k* lies.
#'
#' This function provides a programmatic way to select *k*.
#'  Note that it is highly recommended to look
#' at the graph of BIC versus the numbers of clusters, to understand
#' and validate the programmatic selection. The criteria available in this
#' function are:
#' - "diffNgroup": differences between successive values of the summary
#'  statistics (by default, BIC) are split into two groups using a Ward's
#'  clustering method (see ?hclust), to differentiate sharp decrease from mild
#'  decreases or increases. The retained K is the one before the first group
#'  switch. This criterion appears to work well for island/hierarchical models, and decently
#'  for isolation by distance models, albeit with some unstability. It can be
#'  confounded by an initial, very sharp decrease of the test statistics. IF
#'  UNSURE ABOUT THE CRITERION TO USE, USE THIS ONE.
#' - "min": the model with the minimum summary statistics (as specified by
#'  stat argument, BIC by default) is retained. Is likely to work for simple
#'  island model, using BIC. It is likely to fail in models relating to
#'  stepping stones, where the BIC always decreases (albeit by a small
#'  amount) as K increases. In general, this approach tends to over-estimate
#'  the number of clusters.
#' - "goesup": the selected model is the K after which increasing the number
#' of clusters leads to increasing the summary statistics. Suffers from
#' inaccuracy, since i) a steep decrease might follow a small 'bump' of
#' increase of the statistics, and ii) increase might never happen, or
#' happen after negligible decreases. Is likely to work only for clear-cut
#' island models.
#' - "smoothNgoesup": a variant of "goesup", in which the summary statistics
#' is first smoothed using a lowess approach. Is meant to be more accurate
#' than "goesup" as it is less prone to stopping to small 'bumps' in the
#' decrease of the statistics.
#' - "goodfit": another criterion seeking a good fit with a minimum number
#' of clusters. This approach does not rely on differences between successive
#' statistics, but on absolute fit. It selects the model with the smallest
#' K so that the overall fit is above a given threshold.
#'
#' @references Jombart T, Devillard S and Balloux F (2010) Discriminant analysis of
#' principal components: a new method for the analysis of genetically
#' structured populations. BMC Genetics 11:94. doi:10.1186/1471-2156-11-94
#'
#' @param x a `gt_cluster_pca` object obtained with [gt_cluster_pca()]
#' @param stat a statistics, one of "BIC", "AIC" or "WSS"
#' @param criterion one of "diffNgroup", "min", "goesup", "smoothNgoesup", "goodfit",
#' see details for a discussion of each approach.
#' @param quiet boolean on whether to silence outputting information to the
#'  screen (defaults to FALSE)
#' @returns a 'gt_cluster_pca' object with an added element 'best_k'
#' @export

gt_cluster_pca_best_k <- function(x, stat = c("BIC", "AIC", "WSS"),
                                  criterion = c("diffNgroup", "min", "goesup",
                                                "smoothNgoesup", "goodfit"),
                                  quiet=FALSE){
  if (!inherits(x, "gt_cluster_pca")){
    stop("'x' should be a 'gt_cluster_pcaers' object generated with 'gt_cluster_pca()'")
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
    temp <- stats::cutree(stats::hclust(stats::dist(diff(x$clusters[[stat]])), method="ward.D"), k=2)
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
  x$best_k <- n.clust
  return(x)
}
