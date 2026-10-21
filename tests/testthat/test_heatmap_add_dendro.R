test_that("heatmap_add_dendro adds requested dendrogram panels", {
  toy_mat <- matrix(
    c(
      NA, 1, 2,
      1, NA, 3,
      2, 3, NA
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(
      c("pop_a", "pop_b", "pop_c"),
      c("pop_a", "pop_b", "pop_c")
    )
  )
  class(toy_mat) <- c("pairwise_matrix", class(toy_mat))

  base_plot <- autoplot(toy_mat,
    order = function(x) {
      stats::hclust(stats::as.dist(x))
    }
  ) +
    ggplot2::scale_fill_viridis_c()

  left_plot <- heatmap_add_dendro(
    base_plot,
    side = "left",
    rel_size = 0.2,
    line_color = "firebrick",
    line_size = 1.25
  )

  expect_s3_class(left_plot, "heatmap_dendro")
  expect_true("dend-left" %in% left_plot$gt$layout$name)
  expect_false("dend-top" %in% left_plot$gt$layout$name)

  top_plot <- heatmap_add_dendro(base_plot, side = "top", rel_size = 0.2)

  expect_s3_class(top_plot, "heatmap_dendro")
  expect_true("dend-top" %in% top_plot$gt$layout$name)
  expect_false("dend-left" %in% top_plot$gt$layout$name)

  both_plot <- heatmap_add_dendro(base_plot, side = "both")

  expect_s3_class(both_plot, "heatmap_dendro")
  expect_true(all(c("dend-left", "dend-top") %in% both_plot$gt$layout$name))
})

test_that("heatmap_add_dendro matches the heatmap order with tied distances", {
  # complete linkage has several equally good merges for this matrix
  toy_mat <- matrix(
    c(
      0, 4, 3, 3,
      4, 0, 6, 3,
      3, 6, 0, 3,
      3, 3, 3, 0
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(
      c("pop_a", "pop_b", "pop_c", "pop_d"),
      c("pop_a", "pop_b", "pop_c", "pop_d")
    )
  )
  order_fun <- function(x) {
    stats::hclust(stats::as.dist(x), method = "complete")
  }
  base_plot <- heatmap_pairwise(toy_mat, order = order_fun) +
    ggplot2::scale_fill_viridis_c()
  # the dendrogram should be the clustering used to order the heatmap
  hc <- order_fun(toy_mat)

  # segment coordinates of a dendrogram panel, ignoring drawing order
  seg_coords <- function(grob) {
    segs <- Filter(function(x) inherits(x, "segments"), grob$children)[[1]]
    coords <- sapply(segs[c("x0", "y0", "x1", "y1")], as.numeric)
    coords[do.call(order, as.data.frame(coords)), ]
  }
  dendro_plot <- heatmap_add_dendro(base_plot)
  for (side in c("left", "top")) {
    drawn <- dendro_plot$gt$grobs[[
      which(dendro_plot$gt$layout$name == paste0("dend-", side))
    ]]
    expected <- .dend_grob(hc, side, "grey30", 0.5)
    expect_equal(seg_coords(drawn), seg_coords(expected))
  }

  # plots without a stored clustering fall back to reclustering
  attr(base_plot, "order_obj") <- NULL
  expect_s3_class(heatmap_add_dendro(base_plot), "heatmap_dendro")
})

test_that("heatmap_add_dendro objects can be drawn and saved", {
  toy_mat <- matrix(
    c(
      NA, 1, 2,
      1, NA, 3,
      2, 3, NA
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(
      c("pop_a", "pop_b", "pop_c"),
      c("pop_a", "pop_b", "pop_c")
    )
  )
  class(toy_mat) <- c("pairwise_matrix", class(toy_mat))

  base_plot <- autoplot(toy_mat,
    order = function(x) {
      stats::hclust(stats::as.dist(x))
    }
  ) +
    ggplot2::scale_fill_viridis_c()

  dendro_plot <- heatmap_add_dendro(base_plot)

  expect_invisible(print(dendro_plot))
  expect_no_error(grid::grid.draw(dendro_plot))

  tmp_png <- tempfile(fileext = ".png")

  expect_no_error(
    ggplot2::ggsave(
      filename = tmp_png,
      plot = dendro_plot,
      width = 4,
      height = 4,
      units = "in"
    )
  )
  expect_true(file.exists(tmp_png))
})

test_that("heatmap_add_dendro validates input plots", {
  plain_plot <- heatmap_pairwise(
    matrix(
      c(
        NA, 1, 2,
        1, NA, 3,
        2, 3, NA
      ),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(
        c("pop_a", "pop_b", "pop_c"),
        c("pop_a", "pop_b", "pop_c")
      )
    )
  )

  expect_error(
    heatmap_add_dendro(1),
    "plot must be a ggplot"
  )

  expect_error(
    heatmap_add_dendro(plain_plot),
    "Dendrograms can only be added"
  )
})

test_that("heatmap_add_dendro requires an hclust order function", {
  toy_mat <- matrix(
    c(
      NA, 1, 2,
      1, NA, 3,
      2, 3, NA
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(
      c("pop_a", "pop_b", "pop_c"),
      c("pop_a", "pop_b", "pop_c")
    )
  )
  class(toy_mat) <- c("pairwise_matrix", class(toy_mat))

  base_plot <- autoplot(toy_mat,
    order = function(x) {
      return(list(order = c(2, 1, 3)))
    }
  ) +
    ggplot2::scale_fill_viridis_c()

  expect_error(
    heatmap_add_dendro(
      base_plot
    ),
    "The order function must return an hclust object"
  )
})
