test_that("example_gt functions", {
  test <- example_gt("gen_tbl")
  expect_true(inherits(test, "gen_tbl"))

  test <- example_gt("grouped_gen_tbl")
  expect_true(inherits(test, "gen_tbl"))
  expect_true(inherits(test, "grouped_gen_tbl"))

  test <- example_gt("gen_tbl_sf")
  expect_true(inherits(test, "gen_tbl"))
  expect_true(inherits(test, "sf"))

  test <- example_gt("grouped_gen_tbl_sf")
  expect_true(inherits(test, "gen_tbl"))
  expect_true(inherits(test, "grouped_gen_tbl"))
  expect_true(inherits(test, "sf"))
})
