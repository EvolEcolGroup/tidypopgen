#' Pairwise distance heatmap
#'
#' Plot a square heatmap from either a full symmetric distance matrix or a tidy
#' tibble of distances. Rows and columns can optionally be reordered using
#' either an explicit ordering vector or a function that computes an ordering
#' (e.g. hierarchical clustering). Dendrograms can be added to the heatmap with
#' [heatmap_add_dendro()]. Note that any decoration (e.g. axis labels, fill
#' scale) should be added to the heatmap before adding dendrograms, as the
#' latter will modify the plot structure and make it difficult to add further
#' customization afterwards.
#'
#' @param x A square numeric matrix OR a data.frame/tibble with 3 columns:
#'   names1, names2, value (half-matrix representation).
#'
#' @param order Optional ordering specification. Can be:
#' \itemize{
#'   \item `NULL` (default): preserve the existing order.
#'   \item A character vector containing all row/column names in the desired
#'     order.
#'   \item A function taking the matrix `x` as input and returning either:
#'     \itemize{
#'       \item an ordering vector of indices, or
#'       \item an object with an `order` component (e.g. from
#'         [stats::hclust()]).
#'     }
#' }
#' @seealso [heatmap_add_dendro()]
#' @return A ggplot object.
#'
#' @export
#'
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#'
#' # With a pairwise matrix:
#' fst_dist <- example_gt %>%
#'   group_by(population) %>%
#'   pairwise_pop_fst(method = "Nei87", type = "pairwise")
#'
#' heatmap_pairwise(fst_dist)
#'
#' # Explicit ordering
#' heatmap_pairwise(
#'   fst_dist,
#'   order = rev(rownames(fst_dist))
#' ) +
#' ggplot2::scale_fill_viridis_c()
#'
#' # Hierarchical clustering
#' heatmap_pairwise(
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
#' heatmap_pairwise(fst_tidy)
#'
#' heatmap_pairwise(
#'   fst_tidy,
#'   order = function(x) {
#'     stats::hclust(stats::as.dist(x))
#'   }
#' )
#'
heatmap_pairwise <- function(x,
                             order = NULL) {
  # ============================================================
  # Input validation
  # ============================================================

  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("x must be a matrix or a data frame.")
  }

  # ============================================================
  # CASE 1: long-format input (half matrix, as generated with "tidy")
  # ============================================================

  if (is.data.frame(x)) {
    if (ncol(x) != 3) {
      stop("Data frame must have 3 columns: names_1, names_2, value.")
    }

    if (!grepl("1$", colnames(x)[1]) ||
        !grepl("2$", colnames(x)[2])) {
      stop("First two columns must end in '1' and '2'.")
    }

    df <- x
    colnames(df)[1:3] <- c("from", "to", "value")

    ids <- sort(unique(c(df$from, df$to)))
    n <- length(ids)

    if (nrow(df) != n * (n - 1) / 2) {
      stop("Expected n*(n-1)/2 rows for half-matrix representation.")
    }

    mat <- matrix(
      0,
      n,
      n,
      dimnames = list(ids, ids)
    )

    i <- match(df$from, ids)
    j <- match(df$to, ids)

    mat[cbind(i, j)] <- df$value
    mat[cbind(j, i)] <- df$value

    diag(mat) <- NA

    x <- mat
  }

  # ============================================================
  # CASE 2: matrix input (as generated with "pairwise")
  # ============================================================

  if (is.matrix(x)) {
    if (nrow(x) != ncol(x)) {
      stop("Matrix must be square.")
    }

    if (is.null(rownames(x))) {
      rownames(x) <- seq_len(nrow(x))
    }

    if (is.null(colnames(x))) {
      colnames(x) <- seq_len(ncol(x))
    }

    if (!identical(rownames(x), colnames(x))) {
      stop("Row and column names must match.")
    }
  }

  n <- nrow(x)
  ids <- rownames(x)

  # ============================================================
  # Optional ordering
  # ============================================================

  if (!is.null(order)) {

    # ------------------------------------------------------------
    # Explicit ordering vector
    # ------------------------------------------------------------

    if (is.vector(order) && !is.function(order)) {

      if (length(order) != n) {
        stop("order vector must have length equal to nrow(x).")
      }

      if (!setequal(order, ids)) {
        stop("order vector must contain exactly the row/column names.")
      }

      ord <- match(order, ids)

      # ------------------------------------------------------------
      # Ordering function
      # ------------------------------------------------------------

    } else if (is.function(order)) {

      # make sure there are no NA/NaN values
      d <- stats::as.dist(x)

      if (anyNA(d) || any(is.nan(d))) {
        stop("Cannot compute ordering with NA or NaN values.")
      }

      ord_obj <- order(x)

      # object with $order component
      if (!is.null(ord_obj$order)) {

        ord <- ord_obj$order

        # raw ordering vector
      } else {

        ord <- ord_obj
      }

      if (!is.numeric(ord)) {
        stop("Ordering function must return numeric indices.")
      }

      if (length(ord) != n) {
        stop("Ordering vector has incorrect length.")
      }

      if (!setequal(ord, seq_len(n))) {
        stop("Ordering vector must contain indices 1:n exactly once.")
      }

    } else {

      stop("order must be NULL, a vector, or a function.")
    }

    x <- x[ord, ord]
    ids <- ids[ord]
  }

  # ============================================================
  # Long format for ggplot
  # ============================================================

  df_plot <- expand.grid(
    row = seq_len(n),
    col = seq_len(n)
  )

  df_plot$value <- x[cbind(df_plot$row, df_plot$col)]

  # ============================================================
  # Plot
  # ============================================================
  p <- ggplot2::ggplot(
    df_plot,
    ggplot2::aes(
      x = .data$col,
      y = .data$row,
      fill = .data$value
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::coord_fixed() +
    ggplot2::scale_x_continuous(
      breaks = seq_len(n),
      labels = ids,
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq_len(n),
      labels = ids,
      expand = c(0, 0)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = -60,
        hjust = 0,
        vjust = 1
      )
    )
  attr(p, "order_fun") <- order
  return(p)
}
