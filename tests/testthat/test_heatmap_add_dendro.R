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
