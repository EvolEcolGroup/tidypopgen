test_that("example_gt functions", {
  test <- load_example_gt("gen_tbl")
  expect_true(inherits(test, "gen_tbl"))

  test <- load_example_gt("grouped_gen_tbl")
  expect_true(inherits(test, "gen_tbl"))
  expect_true(inherits(test, "grouped_gen_tbl"))

  test <- load_example_gt("gen_tbl_sf")
  expect_true(inherits(test, "gen_tbl"))
  expect_true(inherits(test, "sf"))

  test <- load_example_gt("grouped_gen_tbl_sf")
  expect_true(inherits(test, "gen_tbl"))
  expect_true(inherits(test, "grouped_gen_tbl"))
  expect_true(inherits(test, "sf"))
})
