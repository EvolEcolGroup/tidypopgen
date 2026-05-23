#' Autoplot for pairwise distance matrices
#'
#' This function provides a convenient way to visualize pairwise distance
#' matrices (e.g., Fst, IBS) as heatmaps using ggplot2. It supports both matrix
#' and tidy tibble formats, and can optionally perform hierarchical clustering
#' for better visualization (no missing values are allowed if clustering).
#'
#' @param object A square numeric matrix OR a data.frame/tibble with 3 columns:
#'  names1, names2, value (half-matrix representation).
#' @param cluster Logical. If TRUE, reorder using hierarchical clustering
#'  (based on [stats::hclust()]).
#' @param hclust_method Character. The agglomeration method for
#'  [stats::hclust()]. It defaults to "complete".
#' @param ... Additional arguments (not currently used). Included for
#' compatibility with autoplot generic.
#' @return A ggplot object.
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#' #'  # With a pairwise matrix:
#' fst_dist <- example_gt %>%
#'   group_by(population) %>%
#'   pairwise_pop_fst(method = "Nei87", type = "pairwise")
#' autoplot(fst_dist, cluster = FALSE)
#' autoplot(fst_dist, cluster = TRUE)
#' #'  # With a tidy tibble of distances:
#' fst_tidy <- example_gt %>%
#'   group_by(population) %>%
#'   pairwise_pop_fst(method = "Nei87", type = "tidy")
#' autoplot(fst_tidy, cluster = FALSE)
#' autoplot(fst_tidy, cluster = TRUE)
#' @name autoplot_pairwise
#' @aliases autoplot_pairwise_matrix autoplot_pairwise_tbl
#' @export
autoplot.pairwise_matrix <- function(object, cluster = FALSE,
                                     hclust_method = "complete",
                                     ...) {
  # check that ellipse is empty
  if (!missing(...)) {
    stop("Additional arguments are not allowed autoplot.pairwise_matrix.")
  }
  heatmap_pairwise(object, cluster = cluster, hclust_method = hclust_method)
}

#' @rdname autoplot_pairwise
#' @export
autoplot.pairwise_tbl <- function(object, cluster = FALSE,
                                  hclust_method = "complete",
                                  ...) {
  if (!missing(...)) {
    stop("Additional arguments are not allowed autoplot.pairwise_matrix.")
  }
  heatmap_pairwise(object, cluster = cluster, hclust_method = hclust_method)
}
