test_that("fct_inorder_base orders levels by first appearance", {
  # character
  f <- fct_inorder_base(c("b", "b", "a", "c", "a"))
  expect_equal(levels(f), c("b", "a", "c"))
  expect_equal(as.character(f), c("b", "b", "a", "c", "a"))
  # factor, including levels that do not appear in the data
  f <- fct_inorder_base(factor(c("b", "a", "b"), levels = c("a", "b", "z")))
  expect_equal(levels(f), c("b", "a", "z"))
  # numeric, which used to become all NA
  f <- fct_inorder_base(c(10, 10, 2, 3, 2))
  expect_equal(levels(f), c("10", "2", "3"))
  expect_equal(as.character(f), c("10", "10", "2", "3", "2"))
  expect_false(anyNA(f))
  # ordered
  expect_true(is.ordered(fct_inorder_base(c(2, 1), ordered = TRUE)))
})
