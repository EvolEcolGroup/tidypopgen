test_that("autoplot works for pairwise matrices and tibbles", {
  # load the example dataset
  example_gt <- load_example_gt("gen_tbl")

  # ============================================================
  # KING
  # ============================================================

  # Matrix form
  king_mat <- pairwise_king(example_gt, as_matrix = TRUE)

  p <- autoplot(king_mat)

  expect_true(inherits(p, "ggplot"))

  # clustering/order function should fail because KING contains NA
  expect_error(
    autoplot(
      king_mat,
      order = function(x) {
        stats::hclust(stats::as.dist(x))
      }
    ),
    "Cannot compute ordering with NA or NaN values."
  )

  # explicit ordering vector
  p <- autoplot(
    king_mat,
    order = rev(rownames(king_mat))
  )

  expect_true(inherits(p, "ggplot"))

  # Tidy tibble form
  king_tbl <- pairwise_king(example_gt, as_matrix = FALSE)

  p <- autoplot(king_tbl)

  expect_true(inherits(p, "ggplot"))

  # break the names of columns
  names(king_tbl)[1] <- "item"

  expect_error(
    autoplot(king_tbl),
    "First two columns must"
  )

  # ============================================================
  # GRM
  # ============================================================

  grm_mat <- example_gt %>%
    pairwise_grm()

  p <- autoplot(grm_mat)

  expect_true(inherits(p, "ggplot"))

  # order using hierarchical clustering
  p <- autoplot(
    grm_mat,
    order = function(x) {
      stats::hclust(
        stats::as.dist(x),
        method = "complete"
      )
    }
  )

  expect_true(inherits(p, "ggplot"))

  # ============================================================
  # IBS
  # ============================================================

  ibs_mat <- example_gt %>%
    pairwise_ibs(as_matrix = TRUE)

  p <- autoplot(ibs_mat)

  expect_true(inherits(p, "ggplot"))

  # explicit ordering
  p <- autoplot(
    ibs_mat,
    order = rownames(ibs_mat)
  )

  expect_true(inherits(p, "ggplot"))

  # invalid ordering vector
  expect_error(
    autoplot(
      ibs_mat,
      order = rownames(ibs_mat)[-1]
    ),
    "order vector must have length"
  )

  # ============================================================
  # Fst
  # ============================================================

  # tidy tibble
  fst_tbl <- example_gt %>%
    dplyr::group_by(population) %>%
    pairwise_pop_fst(method = "Nei87")

  p <- autoplot(fst_tbl)

  expect_true(inherits(p, "ggplot"))

  p <- autoplot(
    fst_tbl,
    order = function(x) {
      stats::hclust(stats::as.dist(x))
    }
  )

  expect_true(inherits(p, "ggplot"))

  # matrix form
  fst_mat <- example_gt %>%
    dplyr::group_by(population) %>%
    pairwise_pop_fst(
      method = "Nei87",
      type = "pairwise"
    )

  p <- autoplot(fst_mat)

  expect_true(inherits(p, "ggplot"))

  p <- autoplot(
    fst_mat,
    order = rev(rownames(fst_mat))
  )

  expect_true(inherits(p, "ggplot"))

  # invalid ordering vector: wrong populations
  expect_error(
    autoplot(
      fst_mat,
      order = c(
        rownames(fst_mat)[-1],
        "not_a_population"
      )
    ),
    "order vector must contain exactly the row/column names"
  )
})

test_that("heatmap_pairwise detects duplicate pairs in tidy input", {
  example_gt <- load_example_gt("gen_tbl")

  fst_tbl <- example_gt %>%
    dplyr::group_by(population) %>%
    pairwise_pop_fst(method = "Nei87")

  # introduce a duplicate pair by rbinding a copy of the first row
  dup_tbl <- rbind(fst_tbl[-nrow(fst_tbl), ], fst_tbl[1, ])

  expect_error(
    heatmap_pairwise(dup_tbl),
    "Duplicate pairs found"
  )
})


test_that("heatmap_pairwise error messages", {
  # load the example dataset
  example_gt <- load_example_gt("gen_tbl")
  # tidy tibble
  fst_tbl <- example_gt %>%
    dplyr::group_by(population) %>%
    pairwise_pop_fst(method = "Nei87")
  # error for additional parameters
  expect_error(
    autoplot(fst_tbl, blah = "blah"),
    "Additional arguments are not allowed"
  )
  # now with a matrix
  # matrix form
  fst_mat <- example_gt %>%
    dplyr::group_by(population) %>%
    pairwise_pop_fst(
      method = "Nei87",
      type = "pairwise"
    )
  expect_error(
    autoplot(fst_mat, blah = "blah"),
    "Additional arguments are not allowed"
  )
  # wrong number of columsn for the tibble
  fst_tbl_2col <- fst_tbl[, 1:2]
  expect_error(
    autoplot(fst_tbl_2col),
    "Data frame must have 3 columns"
  )
  # remove one row
  fst_tbl_missing_row <- fst_tbl[-1, ]
  expect_error(
    autoplot(fst_tbl_missing_row),
    "Expected"
  )
  # create a diagonal entry
  fst_tbl_diag <- fst_tbl
  fst_tbl_diag[1, 1] <- fst_tbl_diag[1, 2]
  expect_error(
    autoplot(fst_tbl_diag),
    "the data.frame should not include values for the diagonal"
  )
})
