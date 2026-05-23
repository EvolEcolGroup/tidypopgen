#' Pairwise distance heatmap
#'
#' Plot a square heatmap from either a full symmetric distance matrix or a tidy
#' tibble of distances. It supports optional hierarchical clustering for
#' row/column reordering using [stats::hclust()].
#'
#' @param x A square numeric matrix OR a data.frame/tibble with 3 columns:
#'   names1, names2, value (half-matrix representation).
#'
#' @param cluster Logical. If TRUE, reorder using hierarchical clustering
#' (based on [stats::hclust()]).
#'
#' @param hclust_method Character. The agglomeration method for
#'   [stats::hclust()]. It defaults to "complete".
#'
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
#' heatmap_pairwise(fst_dist, cluster = FALSE)
#' heatmap_pairwise(fst_dist, cluster = TRUE)
#'
#' # With a tidy tibble of distances:
#' fst_tidy <- example_gt %>%
#'   group_by(population) %>%
#'   pairwise_pop_fst(method = "Nei87", type = "tidy")
#' heatmap_pairwise(fst_tidy, cluster = FALSE)
#' heatmap_pairwise(fst_tidy, cluster = TRUE)
#'
heatmap_pairwise <- function(x,
                             cluster = FALSE,
                             hclust_method = "complete") {
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

    if (!grepl("1$", colnames(x)[1]) || !grepl("2$", colnames(x)[2])) {
      stop("First two columns must end in '1' and '2'.")
    }

    df <- x
    colnames(df)[1:3] <- c("from", "to", "value")

    ids <- sort(unique(c(df$from, df$to)))
    n <- length(ids)

    if (nrow(df) != n * (n - 1) / 2) {
      stop("Expected n*(n-1)/2 rows for half-matrix representation.")
    }

    mat <- matrix(0, n, n,
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
  }

  n <- nrow(x)
  ids <- rownames(x)

  # ============================================================
  # Optional clustering
  # ============================================================

  if (cluster) {
    # make sure there are no na or nan
    if (anyNA(stats::as.dist(x)) || any(is.nan(stats::as.dist(x)))) {
      stop("Cannot cluster with NA or NaN values in the matrix.")
    }
    hc <- stats::hclust(stats::as.dist(x), method = hclust_method)
    ord <- hc$order
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

  ggplot2::ggplot(
    df_plot,
    ggplot2::aes(
      x = .data$col, y = .data$row,
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
      labels = rev(ids),
      expand = c(0, 0)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = -60, hjust = 1)
    )
}
