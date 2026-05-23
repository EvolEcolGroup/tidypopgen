test_that("autoplot works for pairwise matrices and tibbles", {
  # load the example dataset
  example_gt <- load_example_gt("gen_tbl")

  # Compute the KING-robust matrix
  king_mat <- pairwise_king(example_gt, as_matrix = TRUE)
  p <- autoplot(king_mat)
  # this should be a ggplot object
  expect_true(inherits(p, "ggplot"))
  expect_error(
    autoplot(king_mat, cluster = TRUE),
    "Cannot cluster with NA or NaN values in the matrix."
  )
  # Do the same with a tidy tibble
  king_tbl <- pairwise_king(example_gt, as_matrix = FALSE)
  p <- autoplot(king_tbl)
  expect_true(inherits(p, "ggplot"))
  # break the names of columns
  names(king_tbl)[1] <- "item"
  expect_error(
    autoplot(king_tbl),
    "First two columns must"
  )

  # Do the same for GRM
  grm_mat <- example_gt %>% pairwise_grm()
  p <- autoplot(grm_mat)
  expect_true(inherits(p, "ggplot"))

  # IBS
  ibs_mat <- example_gt %>% pairwise_ibs()
  p <- autoplot(ibs_mat)
  expect_true(inherits(p, "ggplot"))

  # Fst
  fst_tbl <- example_gt %>%
    group_by(population) %>%
    pairwise_pop_fst(method = "Nei87")
  p <- autoplot(fst_tbl)
  expect_true(inherits(p, "ggplot"))
  # do the same for a matrix
  fst_mat <- example_gt %>%
    group_by(population) %>%
    pairwise_pop_fst(method = "Nei87", type = "pairwise")
  p <- autoplot(fst_mat)
  expect_true(inherits(p, "ggplot"))
})
