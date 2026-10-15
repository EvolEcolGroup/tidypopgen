#' Autoplot for pairwise distance matrices
#'
#' This function provides a convenient way to visualize pairwise distance
#' matrices (e.g., Fst, IBS) as heatmaps using ggplot2. It supports both matrix
#' and tidy tibble formats, and can optionally reorder rows/columns using either
#' an explicit ordering vector or a function (e.g. hierarchical clustering).
#'
#' @param object A square numeric matrix OR a data.frame/tibble with 3 columns:
#'   names1, names2, value (half-matrix representation).
#'
#' @param order Optional ordering specification. Can be:
#' \itemize{
#'   \item `NULL` (default): preserve the existing order.
#'   \item A character vector containing all row/column names in the desired
#'     order.
#'   \item A function taking the matrix as input and returning either:
#'     \itemize{
#'       \item an ordering vector of indices, or
#'       \item an object with an `order` component (e.g. from
#'         [stats::hclust()]).
#'     }
#' }
#'
#' @param ... Additional arguments (not currently used). Included for
#'   compatibility with autoplot generic.
#'
#' @return A ggplot object.
#'
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # With a pairwise matrix:
#' fst_dist <- example_gt %>%
#'   group_by(population) %>%
#'   pairwise_pop_fst(method = "Nei87", type = "pairwise")
#'
#' autoplot(fst_dist)
#'
#' # Explicit ordering
#' autoplot(
#'   fst_dist,
#'   order = rev(rownames(fst_dist))
#' )
#'
#' # Hierarchical clustering
#' autoplot(
#'   fst_dist,
#'   order = function(x) {
#'     stats::hclust(stats::as.dist(x), method = "complete")
#'   }
#' )
#'
#' # With a tidy tibble of distances:
#' fst_tidy <- example_gt %>%
#'   group_by(population) %>%
#'   pairwise_pop_fst(method = "Nei87", type = "tidy")
#'
#' autoplot(fst_tidy)
#'
#' autoplot(
#'   fst_tidy,
#'   order = function(x) {
#'     stats::hclust(stats::as.dist(x))
#'   }
#' )
#'
#' @name autoplot_pairwise
#' @aliases autoplot_pairwise_matrix autoplot_pairwise_tbl
#' @export
autoplot.pairwise_matrix <- function(object,
                                     order = NULL,
                                     ...) {
  if (!missing(...)) {
    stop(
      "Additional arguments are not allowed in ",
      "autoplot.pairwise_matrix()."
    )
  }

  heatmap_pairwise(object, order = order)
}

#' @rdname autoplot_pairwise
#' @export
autoplot.pairwise_tbl <- function(object,
                                  order = NULL,
                                  ...) {
  if (!missing(...)) {
    stop(
      "Additional arguments are not allowed in ",
      "autoplot.pairwise_tbl()."
    )
  }

  heatmap_pairwise(object, order = order)
}
