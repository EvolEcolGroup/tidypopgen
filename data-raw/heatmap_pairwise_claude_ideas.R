# ggplot2 Heatmap — Modular Implementation
#
# Designed for symmetric distance matrices. A single clustering covers both
# axes, so rows and columns are always reordered identically.
#
# Composable pieces:
#
#   cluster_order(mat)            — hclust object from a distance matrix
#   heatmap_data(mat, hc)         — tidy long data frame (inspectable)
#   gg_heatmap(mat, hc)           — ggplot object; add scales with +
#   scale_fill_heatmap_*()        — named scale_fill wrappers, added with +
#   annotate_dendrogram(plot, hc) — attaches a dendrogram (use |>),
#                                   returns a <heatmap_figure> object
#
# <heatmap_figure> behaves like a ggplot/patchwork object:
#   - bare evaluation and |> pipelines auto-draw (print.heatmap_figure)
#   - ggsave() works directly (grid.draw.heatmap_figure)
#   - chaining annotate_dendrogram() twice works (left then top, or vice-versa)
#
# Alignment note:
#   patchwork aligns bounding boxes, not panels. The heatmap's y-axis labels
#   offset the tile panel, creating a gap and misaligning the top dendrogram.
#   annotate_dendrogram() operates at the gtable level, inserting the dendrogram
#   grob directly into the panel row/column for pixel-perfect alignment.
#
# Typical usage:
#
#   hc <- cluster_order(mat)
#
#   gg_heatmap(mat, hc) |>
#     annotate_dendrogram(hc, "left") |>
#     annotate_dendrogram(hc, "top")
#
#   # Add a scale before piping (+ binds tighter than |>, so parenthesise)
#   (gg_heatmap(mat, hc) + scale_fill_heatmap_sequential(low = "white", high = "darkred")) |>
#     annotate_dendrogram(hc, "left") |>
#     annotate_dendrogram(hc, "top")

library(ggplot2)
library(gtable)
library(grid)
library(scales)


# ══════════════════════════════════════════════════════════════════════════════
# 1.  <heatmap_figure> S3 class
#     A thin wrapper around a gtable that auto-draws on bare evaluation and
#     is accepted by ggsave(), just like a ggplot or patchwork object.
# ══════════════════════════════════════════════════════════════════════════════

new_heatmap_figure <- function(gt) {
  structure(list(gt = gt), class = "heatmap_figure")
}

#' @export
print.heatmap_figure <- function(x, ...) {
  grid.newpage()
  grid.draw(x$gt)
  invisible(x)
}

#' @export
grid.draw.heatmap_figure <- function(x, recording = TRUE) {
  grid.draw(x$gt, recording = recording)
}


# ══════════════════════════════════════════════════════════════════════════════
# 2.  cluster_order() — hierarchical clustering for a distance matrix
# ══════════════════════════════════════════════════════════════════════════════

#' Cluster a symmetric distance matrix
#'
#' Returns an hclust object. Because the input is a distance matrix the same
#' object is used for both axes of the heatmap.
#'
#' @param mat          Numeric symmetric matrix or a dist object.
#' @param dist_method  Distance method passed to dist() (ignored if mat is a dist).
#' @param clust_method Linkage method passed to hclust(). Default "complete".
#' @return An hclust object.
#' @examples
#' hc <- cluster_order(mat)
#' hc <- cluster_order(mat, clust_method = "ward.D2")
cluster_order <- function(mat,
                          dist_method  = "euclidean",
                          clust_method = "complete") {
  d <- if (inherits(mat, "dist")) mat else dist(mat, method = dist_method)
  hclust(d, method = clust_method)
}


# ══════════════════════════════════════════════════════════════════════════════
# 3.  heatmap_data() — matrix → tidy long data frame
# ══════════════════════════════════════════════════════════════════════════════

#' Tidy a matrix into a long data frame ready for ggplot
#'
#' Exported so you can inspect or transform the data before plotting.
#'
#' @param mat  Numeric symmetric matrix.
#' @param hc   hclust object (from cluster_order()). NULL = original order.
#' @return data.frame with columns row (factor), col (factor), value.
heatmap_data <- function(mat, hc = NULL) {
  stopifnot(is.matrix(mat), is.numeric(mat))
  rn <- if (!is.null(rownames(mat))) rownames(mat) else paste0("R", seq_len(nrow(mat)))
  cn <- if (!is.null(colnames(mat))) colnames(mat) else paste0("C", seq_len(ncol(mat)))
  rownames(mat) <- rn
  colnames(mat) <- cn
  if (!is.null(hc)) mat <- mat[hc$order, hc$order, drop = FALSE]
  data.frame(
    row   = factor(rep(rownames(mat), times = ncol(mat)), levels = rev(rownames(mat))),
    col   = factor(rep(colnames(mat), each  = nrow(mat)), levels = colnames(mat)),
    value = as.vector(mat)
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# 4.  gg_heatmap() — returns a plain ggplot
# ══════════════════════════════════════════════════════════════════════════════

#' Create a ggplot2 heatmap
#'
#' Returns a plain ggplot. Add scales with + and dendrogram panels with
#' |> annotate_dendrogram(), which returns a <heatmap_figure> that prints and
#' saves like a ggplot.
#'
#' @param mat          Numeric symmetric matrix.
#' @param hc           hclust object (from cluster_order()). NULL = original order.
#' @param na_color     Fill colour for NA cells.
#' @param cell_color   Colour of cell borders (NA = none).
#' @param show_values  Overlay cell values as text?
#' @param value_format sprintf format string for cell labels.
#' @return A ggplot object.
#' @examples
#' hc <- cluster_order(mat)
#'
#' gg_heatmap(mat, hc)
#'
#' gg_heatmap(mat, hc) + scale_fill_heatmap_sequential(low = "white", high = "darkred")
#'
#' gg_heatmap(mat, hc) |>
#'   annotate_dendrogram(hc, "left") |>
#'   annotate_dendrogram(hc, "top")
gg_heatmap <- function(mat,
                       hc           = NULL,
                       na_color     = "grey85",
                       cell_color   = NA,
                       show_values  = FALSE,
                       value_format = "%.2f") {
  df <- heatmap_data(mat, hc)
  
  p <- ggplot(df, aes(x = col, y = row, fill = value)) +
    geom_tile(color = cell_color, linewidth = 0.4) +
    scale_fill_gradient2(
      low      = "#2166AC",
      mid      = "#F7F7F7",
      high     = "#D6604D",
      midpoint = mean(mat, na.rm = TRUE),
      na.value = na_color,
      name     = "Value"
    ) +
    coord_fixed() +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      axis.text.x       = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid        = element_blank(),
      legend.key.height = unit(1.2, "cm")
    )
  
  if (show_values) {
    rng     <- range(df$value, na.rm = TRUE)
    norm    <- rescale(df$value, from = rng)
    df$tcol <- ifelse(is.na(norm), "black", ifelse(norm > 0.55, "white", "black"))
    p <- p + geom_text(
      data    = df,
      mapping = aes(label = ifelse(is.na(value), "NA", sprintf(value_format, value))),
      colour  = df$tcol,
      size    = 3
    )
  }
  
  p
}


# ══════════════════════════════════════════════════════════════════════════════
# 5.  Scale wrappers
# ══════════════════════════════════════════════════════════════════════════════

#' Diverging blue–white–red fill scale
#' @param low,mid,high Endpoint and midpoint colours.
#' @param midpoint     Value mapped to mid colour (default 0).
#' @param na.value     Colour for NA cells.
#' @param ...          Forwarded to scale_fill_gradient2().
scale_fill_heatmap_diverging <- function(low      = "#2166AC",
                                         mid      = "#F7F7F7",
                                         high     = "#D6604D",
                                         midpoint = 0,
                                         na.value = "grey85",
                                         ...) {
  scale_fill_gradient2(low = low, mid = mid, high = high,
                       midpoint = midpoint, na.value = na.value, ...)
}

#' Sequential single-hue fill scale
#' @param low,high Colour endpoints.
#' @param na.value Colour for NA cells.
#' @param ...      Forwarded to scale_fill_gradient().
scale_fill_heatmap_sequential <- function(low      = "#EFF3FF",
                                          high     = "#084594",
                                          na.value = "grey85",
                                          ...) {
  scale_fill_gradient(low = low, high = high, na.value = na.value, ...)
}

#' Viridis fill scale
#' @param option Viridis palette: "viridis", "magma", "plasma", "inferno", "cividis".
#' @param ...    Forwarded to scale_fill_viridis_c().
scale_fill_heatmap_viridis <- function(option = "viridis", ...) {
  scale_fill_viridis_c(option = option, ...)
}


# ══════════════════════════════════════════════════════════════════════════════
# 6.  annotate_dendrogram() — inserts a dendrogram panel via gtable
# ══════════════════════════════════════════════════════════════════════════════

#' Attach a dendrogram panel to a heatmap
#'
#' Designed for use with |>. Accepts a ggplot or <heatmap_figure> and returns
#' a <heatmap_figure> with the dendrogram panel inserted directly adjacent to
#' the heatmap panel — pixel-perfect alignment regardless of axis label width.
#' Can be chained to add dendrograms on both sides.
#'
#' @param plot       A ggplot (from gg_heatmap()) or <heatmap_figure> from a
#'                   prior annotate_dendrogram() call.
#' @param hc         hclust object (from cluster_order()).
#' @param side       "left" (row dendrogram) or "top" (column dendrogram).
#' @param rel_size   Size of the dendrogram relative to the panel (null unit).
#' @param line_color Segment colour.
#' @param line_size  Segment linewidth.
#' @return A <heatmap_figure> object.
#' @examples
#' hc <- cluster_order(mat)
#'
#' gg_heatmap(mat, hc) |> annotate_dendrogram(hc, "left")
#'
#' gg_heatmap(mat, hc) |>
#'   annotate_dendrogram(hc, "left") |>
#'   annotate_dendrogram(hc, "top")
annotate_dendrogram <- function(plot,
                                hc,
                                side       = c("left", "top"),
                                rel_size   = 0.15,
                                line_color = "grey30",
                                line_size  = 0.5) {
  side <- match.arg(side)
  gt   <- if (inherits(plot, "ggplot"))         ggplotGrob(plot) else
    if (inherits(plot, "heatmap_figure"))  plot$gt          else
      stop("plot must be a ggplot or heatmap_figure")
  pp   <- .panel_pos(gt)
  grob <- .dend_grob(hc, side, line_color, line_size)
  
  if (side == "left") {
    gt <- gtable_add_cols(gt, unit(rel_size, "null"), pos = pp$l - 1)
    gt <- gtable_add_grob(gt, grob,
                          t = pp$t, l = pp$l, b = pp$b, r = pp$l,
                          name = "dend-left")
  } else {
    gt <- gtable_add_rows(gt, unit(rel_size, "null"), pos = pp$t - 1)
    gt <- gtable_add_grob(gt, grob,
                          t = pp$t, l = pp$l, b = pp$t, r = pp$r,
                          name = "dend-top")
  }
  
  new_heatmap_figure(gt)
}


# ══════════════════════════════════════════════════════════════════════════════
# 7.  Internal helpers
# ══════════════════════════════════════════════════════════════════════════════

# Row/col indices of the panel cell in a gtable.
.panel_pos <- function(gt) {
  gt$layout[gt$layout$name == "panel", c("t", "l", "b", "r")]
}

# Build a dendrogram ggplot and return just its panel grob.
# Each leaf is mapped to its display rank in hc$order, matching the heatmap
# factor axis exactly.
.dend_grob <- function(hc, side, line_color, line_size) {
  seg_df <- .hclust_segments(hc)
  n      <- length(hc$order)
  
  if (side == "left") {
    p <- ggplot(seg_df, aes(x = -y, xend = -yend, y = x, yend = xend)) +
      geom_segment(color = line_color, linewidth = line_size) +
      scale_y_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
      scale_x_continuous(expand = expansion(mult = c(0.02, 0.1))) +
      theme_void()
  } else {
    p <- ggplot(seg_df, aes(x = x, xend = xend, y = y, yend = yend)) +
      geom_segment(color = line_color, linewidth = line_size) +
      scale_x_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
      theme_void()
  }
  
  gt <- ggplotGrob(p)
  gt$grobs[[which(gt$layout$name == "panel")]]
}

# Convert an hclust object to a segment data frame.
# Leaf k is placed at position rank(k, hc$order): position 1 = first
# displayed label, n = last — matching the heatmap factor axis exactly.
.hclust_segments <- function(hc) {
  n      <- length(hc$order)
  n_node <- nrow(hc$merge)
  height <- hc$height
  
  leaf_pos           <- integer(n)
  leaf_pos[hc$order] <- seq_len(n)
  
  pos <- numeric(n + n_node)
  pos[seq_len(n)] <- leaf_pos
  
  segs <- vector("list", n_node)
  for (i in seq_len(n_node)) {
    l  <- hc$merge[i, 1];  r  <- hc$merge[i, 2]
    lx <- if (l < 0) pos[-l] else pos[n + l]
    rx <- if (r < 0) pos[-r] else pos[n + r]
    lh <- if (l < 0) 0       else height[l]
    rh <- if (r < 0) 0       else height[r]
    ch <- height[i]
    pos[n + i] <- (lx + rx) / 2
    segs[[i]]  <- data.frame(
      x    = c(lx, lx, rx), xend = c(lx, rx, rx),
      y    = c(lh, ch, rh), yend = c(ch, ch, ch)
    )
  }
  do.call(rbind, segs)
}


# ══════════════════════════════════════════════════════════════════════════════
# USAGE EXAMPLES
# ══════════════════════════════════════════════════════════════════════════════
if (FALSE) {
  
  set.seed(42)
  raw <- matrix(rnorm(100), nrow = 10)
  mat <- as.matrix(dist(raw))
  dimnames(mat) <- list(paste0("Item", 1:10), paste0("Item", 1:10))
  
  hc <- cluster_order(mat)
  
  # 1. Bare heatmap — ggplot, prints automatically
  gg_heatmap(mat, hc)
  
  # 2. Custom scale
  gg_heatmap(mat, hc) +
    scale_fill_heatmap_sequential(low = "white", high = "darkred", name = "Distance")
  
  # 3. Viridis scale
  gg_heatmap(mat, hc) +
    scale_fill_heatmap_viridis(option = "magma", name = "Distance")
  
  # 4. Left dendrogram — heatmap_figure, prints automatically
  gg_heatmap(mat, hc) |>
    annotate_dendrogram(hc, "left")
  
  # 5. Both dendrograms
  gg_heatmap(mat, hc) |>
    annotate_dendrogram(hc, "left") |>
    annotate_dendrogram(hc, "top")
  
  # 6. Scale + both dendrograms (parentheses needed: + binds tighter than |>)
  (gg_heatmap(mat, hc) +
      scale_fill_heatmap_sequential(low = "white", high = "darkred")) |>
    annotate_dendrogram(hc, "left") |>
    annotate_dendrogram(hc, "top")
  
  # 7. Cell values + ward linkage + styled dendrogram
  hc_ward <- cluster_order(mat, clust_method = "ward.D2")
  (gg_heatmap(mat, hc_ward,
              show_values  = TRUE,
              value_format = "%.1f",
              cell_color   = "white") +
      scale_fill_heatmap_sequential(low = "white", high = "darkred")) |>
    annotate_dendrogram(hc_ward, "left", rel_size = 0.12, line_color = "steelblue")
  
  # 8. ggsave — heatmap_figure accepted directly
  fig <- (gg_heatmap(mat, hc) +
            scale_fill_heatmap_sequential(low = "white", high = "darkred")) |>
    annotate_dendrogram(hc, "left") |>
    annotate_dendrogram(hc, "top")
  ggsave("heatmap.pdf", fig, width = 8, height = 7)
}