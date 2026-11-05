test_that("annotate_group_info ignores unused group levels", {
  q_tbl <- tibble::tibble(
    id = rep(c("a", "b", "c", "d"), each = 2),
    group = factor(
      rep(c("P2", "P2", "P1", "P1"), each = 2),
      levels = c("P1", "P2", "P3")
    ),
    q = rep(c("Q1", "Q2"), 4),
    percentage = 0.5
  )
  plt <- ggplot2::ggplot(
    q_tbl,
    ggplot2::aes(x = id, y = percentage, fill = q)
  ) +
    ggplot2::geom_col()
  annotation <- annotate_group_info(q_tbl, plt)
  # one separator at the end of each group, and no phantom group P3
  expect_equal(annotation[[1]]$data$x, c(2.5, 4.5))
  expect_equal(annotation[[2]]$breaks, c("a", "c"))
  expect_equal(annotation[[2]]$labels, c("P2", "P1"))
})

test_that("annotate_group_info labels groups of a single individual", {
  # groups of one individual, first and in the middle
  q_tbl <- tibble::tibble(
    id = rep(c("a", "b", "c", "d", "e", "f"), each = 2),
    group = rep(c("P1", "P2", "P2", "P3", "P4", "P4"), each = 2),
    q = rep(c("Q1", "Q2"), 6),
    percentage = 0.5
  )
  plt <- ggplot2::ggplot(
    q_tbl,
    ggplot2::aes(x = id, y = percentage, fill = q)
  ) +
    ggplot2::geom_col()
  annotation <- annotate_group_info(q_tbl, plt)
  expect_equal(annotation[[1]]$data$x, c(1.5, 3.5, 4.5, 6.5))
  # each label is on an individual of its group (the left one of the middle
  # two for an even number of individuals, as above)
  expect_equal(annotation[[2]]$breaks, c("a", "b", "d", "e"))
  expect_equal(annotation[[2]]$labels, c("P1", "P2", "P3", "P4"))
  # and the plot can be built
  expect_no_error(ggplot2::ggplot_build(plt + annotation))
})
