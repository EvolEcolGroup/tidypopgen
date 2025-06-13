#' Run K-clustering on principal components
#'
#' This function implements the clustering procedure used in Discriminant
#' Analysis of Principal Components (DAPC, Jombart et al. 2010). This procedure
#' consists in running successive K-means with an increasing number of clusters
#' (k), after transforming data using a principal component analysis (PCA). For
#' each model, several statistical measures of goodness of fit are computed,
#' which allows to choose the optimal k using the function
#' [gt_cluster_pca_best_k()]. See details for a description of how to select the
#' optimal k and vignette("adegenet-dapc") for a tutorial.
#' @param  x a `gt_pca` object returned by one of the `gt_pca_*` functions.
#' @param n_pca number of principal components to be fed to the LDA.
#' @param k_clusters number of clusters to explore, either a single value, or a
#'   vector of length 2 giving the minimum and maximum (e.g. 1:5). If left NULL,
#'   it will use 1 to the number of pca components divided by 10 (a reasonable
#'   guess).
#' @param method either 'kmeans' or 'ward'
#' @param n_iter number of iterations for kmeans (only used if
#'   `method="kmeans"`)
#' @param n_start number of starting points for kmeans (only used if
#'   `method="kmeans"`)
#' @param quiet boolean on whether to silence outputting information to the
#'   screen (defaults to FALSE)
#' @returns a `gt_cluster_pca` object, which is a subclass of `gt_pca` with an
#'   additional element 'cluster', a list with elements:
#' - 'method' the clustering method (either kmeans or ward)
#' - 'n_pca' number of principal components used for clustering
#' - 'k' the k values explored by the function
#' - 'WSS' within sum of squares for each k
#' - 'AIC' the AIC for each k
#' - 'BIC' the BIC for each k
#' - 'groups' a list, with each element giving the group assignments
#'    for a given k
#' @export
#' @examples
#' # Create a gen_tibble of lobster genotypes
#' bed_file <-
#'   system.file("extdata", "lobster", "lobster.bed", package = "tidypopgen")
#' lobsters <- gen_tibble(bed_file,
#'   backingfile = tempfile("lobsters"),
#'   quiet = TRUE
#' )
#'
#' # Remove monomorphic loci and impute
#' lobsters <- lobsters %>% select_loci_if(loci_maf(genotypes) > 0)
#' lobsters <- gt_impute_simple(lobsters, method = "mode")
#'
#' # Create PCA object
#' pca <- gt_pca_partialSVD(lobsters)
#'
#' # Run clustering on the first 10 PCs
#' gt_cluster_pca(
#'   x = pca,
#'   n_pca = 10,
#'   k_clusters = c(1, 5),
#'   method = "kmeans",
#'   n_iter = 1e5,
#'   n_start = 10,
#'   quiet = FALSE
#' )
#'
#' # Alternatively, use method "ward"
#' gt_cluster_pca(
#'   x = pca,
#'   n_pca = 10,
#'   k_clusters = c(1, 5),
#'   method = "ward",
#'   quiet = FALSE
#' )
#'
gt_cluster_pca <- function(
    x = NULL,
    n_pca = NULL,
    k_clusters = c(1, round(nrow(x$u) / 10)),
    method = c("kmeans", "ward"),
    n_iter = 1e5,
    n_start = 10,
    quiet = FALSE) {
  if (is.null(x) || !inherits(x, "gt_pca")) {
    stop("'x' should be a 'gt_pca' object")
  }

  ## CHECKS ##
  method <- match.arg(method)

  ## SOME GENERAL VARIABLES ##
  n <- nrow(x$u)
  # if we don't give the number of pca, then use them all
  if (is.null(n_pca)) {
    n_pca <- length(x$d)
    if (!quiet) {
      message("'n_pca' was not set: all ", n_pca, " components were used")
    }
  }

  # unpack number of clusters
  if (length(k_clusters) == 1) {
    nb_clust <- k_clusters
  } else if (length(k_clusters) == 2) {
    nb_clust <- k_clusters[1]:k_clusters[2]
  } else {
    stop(paste(
      "'k_clusters' should be either a single value, or the",
      "minimum and maximum to be tested"
    ))
  }

  # set up list to store all results
  cluster_list <- list(
    method = method,
    n_pca = n_pca,
    k = nb_clust,
    WSS = numeric(0),
    AIC = numeric(0),
    BIC = numeric(0),
    groups = list()
  )

  # cluster of 1 is a special case
  if (nb_clust[1] == 1) {
    nb_clust <- nb_clust[-1]
    add_clust_of_1 <- TRUE
  } else {
    add_clust_of_1 <- FALSE
  }
  # vector to store within sum of squares
  wss <- numeric(0)
  # get the scores
  x_scores <- sweep(x$u, 2, x$d, "*")[, seq_len(n_pca)]

  # TODO this is an obvious place to parallelise
  for (i in seq_along(nb_clust)) {
    if (method == "kmeans") {
      ## kmeans clustering (original method)
      groups_assignments <- stats::kmeans(
        x_scores,
        centers = nb_clust[i],
        iter.max = n_iter,
        nstart = n_start
      )$cluster
    } else {
      ## ward clustering
      groups_assignments <- stats::cutree(
        stats::hclust(stats::dist(x_scores)^2, method = "ward.D2"),
        k = nb_clust[i]
      ) # nolint
    }
    wss[i] <- compute_wss(x_scores, groups_assignments)
    cluster_list$groups[[i]] <- groups_assignments
  }
  # compute statistics of goodnees of fit
  if (add_clust_of_1) {
    wss_ori <- sum(apply(x_scores, 2, function(v) sum((v - mean(v))^2)))
    wss <- c(wss_ori, wss)
    # add the classification for 1 cluster (they are all 1!)
    cluster_list$groups <- append(
      cluster_list$groups,
      list(stats::setNames(rep(1, n), names(cluster_list$groups[[2]]))),
      after = 0
    )
    nb_clust <- c(1, nb_clust)
  }
  cluster_list$AIC <- n * log(wss / n) + 2 * nb_clust
  cluster_list$BIC <- n * log(wss / n) + log(n) * nb_clust
  cluster_list$WSS <- wss

  names(cluster_list$groups) <- nb_clust

  x$clusters <- cluster_list
  class(x) <- c("gt_cluster_pca", class(x))
  return(x)
}

## Compute within sum of squares from a matrix 'x' and a factor 'f'
compute_wss <- function(x, f) {
  x_group_mean <- apply(x, 2, tapply, f, mean)
  sum((x - x_group_mean[as.character(f), ])^2)
}


#' Autoplots for `gt_cluster_pca` objects
#'
#' For `gt_cluster_pca`, `autoplot` produces a plot of a metric of choice
#' ('BIC', 'AIC' or 'WSS') against the number of clusters (*k*). This plot is
#' can be used to infer the best value of *k*, which corresponds to the smallest
#' value of the metric (the minimum in an 'elbow' shaped curve). In some cases,
#' there is not 'elbow' and the metric keeps decreasing with increasing *k*; in
#' such cases, it is customary to choose the value of *k* at which the decrease
#' in the metric reaches as plateau. For a programmatic way of choosing
#' *k*, use [gt_cluster_pca_best_k()].
#'
#' `autoplot` produces simple plots to quickly inspect an object. They are not
#' customisable; we recommend that you use `ggplot2` to produce publication
#' ready plots.
#'
#' @param object an object of class `gt_dapc`
#' @param metric the metric to plot on the y axis, one of 'BIC', 'AIC', or
#'   'WSS' (with sum of squares)
#' @param ... not currently used.
#' @returns a `ggplot2` object
#' @rdname autoplot_gt_cluster_pca
#' @export
#' @examples
#' # Create a gen_tibble of lobster genotypes
#' bed_file <-
#'   system.file("extdata", "lobster", "lobster.bed", package = "tidypopgen")
#' lobsters <- gen_tibble(bed_file,
#'   backingfile = tempfile("lobsters"),
#'   quiet = TRUE
#' )
#'
#' # Remove monomorphic loci and impute
#' lobsters <- lobsters %>% select_loci_if(loci_maf(genotypes) > 0)
#' lobsters <- gt_impute_simple(lobsters, method = "mode")
#'
#' # Create PCA object
#' pca <- gt_pca_partialSVD(lobsters)
#'
#' # Run clustering on the first 10 PCs
#' cluster_pca <- gt_cluster_pca(
#'   x = pca,
#'   n_pca = 10,
#'   k_clusters = c(1, 5),
#'   method = "kmeans",
#'   n_iter = 1e5,
#'   n_start = 10,
#'   quiet = FALSE
#' )
#'
#' # Autoplot BIC
#' autoplot(cluster_pca, metric = "BIC")
#'
#' # # Autoplot AIC
#' autoplot(cluster_pca, metric = "AIC")
#'
#' # # Autoplot WSS
#' autoplot(cluster_pca, metric = "WSS")
autoplot.gt_cluster_pca <- function(
    object,
    metric = c("BIC", "AIC", "WSS"),
    ...) {
  metric <- match.arg(metric)
  # create a small tibble with the data of interest
  metric_tbl <- tibble::tibble(
    metric = object$clusters[[metric]],
    k = object$clusters$k
  )
  ggplot2::ggplot(
    data = metric_tbl,
    ggplot2::aes(x = .data$k, y = .data$metric)
  ) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::labs(y = metric)
}


# a print method
#' @method print gt_cluster_pca
#' @export
print.gt_cluster_pca <- function(x, ...) {
  cat(" === Clustering of PCs of gen_tibble object ===")
  cat(
    "\nMethod ($clusters$method): ",
    x$clusters$method
  )
  cat("\n")
  cat(
    "\nN of PCs ($clusters$n_pca): ",
    x$clusters$n_pca
  )
  cat("\n")
  cat(
    "\nK for clustering ($clusters$k):",
    min(x$clusters$k), max(x$clusters$k)
  )
  cat("\n")
  cat(
    "\nThe clustering information is in the slot $clusters;"
  )
  cat("\nother slots are the same as in a gt_pca object used for clustering\n")
  cat("\n")
  invisible(x)
}
