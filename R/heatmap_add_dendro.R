#' Add dendrograms to a heatmap of pairwise distances
#'
#' This function adds dendrograms to a heatmap of pairwise distances, generated
#' with `heatmap_pairwise()`. The heatmap must have been created with an `order`
#' function based on `hclust()`, and the dendrograms will be based on the same
#' clustering. Dendrograms can be added to the left (row) and/or top (column) of
#' the heatmap. Note that the heatmap needs to be formatted first, as it can not
#' be modified easily once the dendrogram panels have been added.
#'
#' To correctly align the dendrogram with the labels of the heatmap, we use a
#' `gtable` to insert the dendrogram panel directly adjacent to the heatmap
#' panel. This means that the dendrogram will be perfectly aligned with the
#' heatmap cells regardless of axis label width. However, manipulating
#' individual plots in a gtable is not easy, so the heatmap should be formatted
#' fully before the dendrograms are added.
#'
#' @param plot A ggplot from `heatmap_pairwise()`, with an order function based
#'   on `hclust()`.
#' @param side "left" (row dendrogram), "top" (column dendrogram), or "both".
#' @param rel_size Size of the dendrogram relative to the panel (null unit).
#' @param line_color Segment colour.
#' @param line_size Segment linewidth.
#' @return A `heatmap_dendro` object (a subclass of `gtable`), which auto-draws
#'   on `print()` and is accepted by `plot()` and `ggsave()`, just like a
#'   `ggplot` or `patchwork` object.
#' @export
#' @examples
#' example_gt <- load_example_gt("gen_tbl")
#' # With a pairwise matrix:
#' fst_dist <- example_gt %>%
#'   group_by(population) %>%
#'   pairwise_pop_fst(method = "Nei87", type = "pairwise")
#' mat_plot <- heatmap_pairwise(fst_dist,
#'   order = function(x) {
#'     stats::hclust(stats::as.dist(x))
#'   }) +
#' ggplot2::scale_fill_viridis_c()
#' heatmap_add_dendro(mat_plot)

heatmap_add_dendro <- function(plot,
                               side = c("both", "left", "top"),
                               rel_size = 0.15,
                               line_color = "grey30",
                               line_size = 0.5) {
  # check that plot is a ggplot
  if (!inherits(plot, "ggplot")) {
    stop("plot must be a ggplot or heatmap_dendro")
  }
  # check that the plot has an "order_fun" attribute, and that it is a function
  if (!is.function(attr(plot, "order_fun"))) {
    stop(
      "Dendrograms can only be added to a plot generated with ",
      "`heatmap_pairwise()` and an order function based on hclust ."
    )
  }

  side <- match.arg(side)

  #############################################
  # get data
  #############################################
  df <- plot@data
  # determine matrix size
  n <- max(c(df$row, df$col))
  # initialise matrix
  mat <- matrix(
    NA_real_,
    nrow = n,
    ncol = n
  )
  # fill matrix
  mat[cbind(df$row, df$col)] <- df$value

  #############################################
  # cluster the data
  #############################################
  cluster_fun <- attr(plot, "order_fun")
  hc <- cluster_fun(mat)

  # validate that the clustering function returned an hclust object
  if (!inherits(hc, "hclust")) {
    stop(
      "The order function must return an hclust object for ",
      "heatmap_add_dendro()."
    )
  }

  if (side == "left" || side == "both") {
    plot <- annotate_dendrogram(
      plot, hc, "left", rel_size, line_color, line_size
    )
  }
  if (side == "top" || side == "both") {
    plot <- annotate_dendrogram(
      plot, hc, "top", rel_size, line_color, line_size
    )
  }
  return(plot)
}


#############################################
# <heatmap_dendro> S3 class
#     A thin wrapper around a gtable that auto-draws on bare evaluation and
#     is accepted by ggsave(), just like a ggplot or patchwork object.
#############################################

new_heatmap_dendro <- function(gt) {
  structure(list(gt = gt), class = c("heatmap_dendro", "gtable"))
}

#' @export
print.heatmap_dendro <- function(x, ...) {
  grid::grid.newpage()
  grid::grid.draw(x$gt)
  invisible(x)
}

#' @importFrom grid grid.draw
#' @export
grid.draw.heatmap_dendro <- function(x, recording = TRUE) {
  grid::grid.draw(x$gt, recording = recording)
}


#############################################
# annotate_dendrogram() — inserts a dendrogram panel via gtable
#############################################

#' Attach a dendrogram panel to a heatmap
#'
#' Designed for use with |>. Accepts a ggplot or <heatmap_dendro> and returns
#' a <heatmap_dendro> with the dendrogram panel inserted directly adjacent to
#' the heatmap panel — pixel-perfect alignment regardless of axis label width.
#' Can be chained to add dendrograms on both sides.
#'
#' @param plot A ggplot (from heatmap_pairwise()) or <heatmap_dendro> from a
#'   prior annotate_dendrogram() call.
#' @param hc hclust object (from order functions that return hclust).
#' @param side "left" (row dendrogram) or "top" (column dendrogram).
#' @param rel_size Size of the dendrogram relative to the panel (null unit).
#' @param line_color Segment colour.
#' @param line_size Segment linewidth.
#' @return A <heatmap_dendro> object.
#' @keywords internal
annotate_dendrogram <- function(plot,
                                hc,
                                side = c("left", "top"),
                                rel_size = 0.15,
                                line_color = "grey30",
                                line_size = 0.5) {
  side <- match.arg(side)
  if (inherits(plot, "ggplot")) {
    gt <- ggplot2::ggplotGrob(plot)
  } else if (inherits(plot, "heatmap_dendro")) {
    gt <- plot$gt
  } else {
    stop("plot must be a ggplot or heatmap_dendro")
  }
  pp <- .panel_pos(gt)
  grob <- .dend_grob(hc, side, line_color, line_size)

  if (side == "left") {
    gt <- gtable::gtable_add_cols(
      gt, ggplot2::unit(rel_size, "null"), pos = pp$l - 1
    )
    gt <- gtable::gtable_add_grob(
      gt, grob,
      t = pp$t, l = pp$l, b = pp$b, r = pp$l,
      name = "dend-left"
    )
  } else {
    gt <- gtable::gtable_add_rows(
      gt, ggplot2::unit(rel_size, "null"), pos = pp$t - 1
    )
    gt <- gtable::gtable_add_grob(
      gt, grob,
      t = pp$t, l = pp$l, b = pp$t, r = pp$r,
      name = "dend-top"
    )
  }

  new_heatmap_dendro(gt)
}


#############################################
# Internal helpers
#############################################

# Row/col indices of the panel cell in a gtable.
.panel_pos <- function(gt) {
  gt$layout[gt$layout$name == "panel", c("t", "l", "b", "r")]
}

# Build a dendrogram ggplot and return just its panel grob.
# Each leaf is mapped to its display rank in hc$order, matching the heatmap
# factor axis exactly.
.dend_grob <- function(hc, side, line_color, line_size) {
  seg_df <- .hclust_segments(hc)
  n <- length(hc$order)

  if (side == "left") {
    p <- ggplot2::ggplot(
      seg_df,
      ggplot2::aes(
        x = -.data$y,
        xend = -.data$yend,
        y = .data$x,
        yend = .data$xend
      )
    ) +
      ggplot2::geom_segment(color = line_color, linewidth = line_size) +
      ggplot2::scale_y_continuous(
        limits = c(0.5, n + 0.5), expand = c(0, 0)
      ) +
      ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(mult = c(0.02, 0.1))
      ) +
      ggplot2::theme_void()
  } else {
    p <- ggplot2::ggplot(
      seg_df,
      ggplot2::aes(
        x = .data$x,
        xend = .data$xend,
        y = .data$y,
        yend = .data$yend
      )
    ) +
      ggplot2::geom_segment(color = line_color, linewidth = line_size) +
      ggplot2::scale_x_continuous(
        limits = c(0.5, n + 0.5), expand = c(0, 0)
      ) +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = c(0, 0.1))
      ) +
      ggplot2::theme_void()
  }

  gt <- ggplot2::ggplotGrob(p)
  gt$grobs[[which(gt$layout$name == "panel")]]
}

# Convert an hclust object to a segment data frame.
# Leaf k is placed at position rank(k, hc$order): position 1 = first
# displayed label, n = last — matching the heatmap factor axis exactly.
.hclust_segments <- function(hc) {
  n <- length(hc$order)
  n_node <- nrow(hc$merge)
  height <- hc$height

  leaf_pos <- integer(n)
  leaf_pos[hc$order] <- seq_len(n)

  pos <- numeric(n + n_node)
  pos[seq_len(n)] <- leaf_pos

  segs <- vector("list", n_node)
  for (i in seq_len(n_node)) {
    l <- hc$merge[i, 1]
    r <- hc$merge[i, 2]
    lx <- if (l < 0) pos[-l] else pos[n + l]
    rx <- if (r < 0) pos[-r] else pos[n + r]
    lh <- if (l < 0) 0 else height[l]
    rh <- if (r < 0) 0 else height[r]
    ch <- height[i]
    pos[n + i] <- (lx + rx) / 2
    segs[[i]] <- data.frame(
      x = c(lx, lx, rx), xend = c(lx, rx, rx),
      y = c(lh, ch, rh), yend = c(ch, ch, ch)
    )
  }
  do.call(rbind, segs)
}
