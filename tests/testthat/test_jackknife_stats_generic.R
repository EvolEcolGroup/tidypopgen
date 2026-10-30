test_that("jackknife_stats_generic works with mean, sum and custom", {
  test_gt <- load_example_gt("gen_tbl_gen_dist")
  loci <- show_loci(test_gt)
  x <- c(0.1, 0.4, NA, 0.3, 0.2, 0.6, 0.5, 0.1)
  jack_mean <- jackknife_stats_generic(x,
    loci_table = loci, operator = "mean", block_size = 2
  )
  jack_custom <- jackknife_stats_generic(x,
    loci_table = loci, operator = "custom", block_size = 2,
    f = mean, na.rm = TRUE
  )
  expect_equal(jack_mean, jack_custom)
  expect_equal(jack_mean$estimate, mean(x, na.rm = TRUE))
  # the second block (a single NA locus) is dropped
  expect_equal(jack_mean$n_blocks, 4L)

  jack_sum <- jackknife_stats_generic(x,
    loci_table = loci, operator = "sum", block_size = 1, size_unit = "cM"
  )
  jack_sum_custom <- jackknife_stats_generic(x,
    loci_table = loci, operator = "custom", block_size = 1,
    size_unit = "cM", f = sum, na.rm = TRUE
  )
  expect_equal(jack_sum, jack_sum_custom)
  expect_equal(jack_sum$estimate, sum(x, na.rm = TRUE))
})

test_that("jackknife_stats_generic checks its inputs", {
  test_gt <- load_example_gt("gen_tbl_gen_dist")
  loci <- show_loci(test_gt)
  expect_error(
    jackknife_stats_generic(1:3, loci_table = loci, block_size = 2),
    "loci_table must have the same number of rows as x"
  )
  expect_error(
    jackknife_stats_generic(1:8,
      loci_table = loci, operator = "custom",
      block_size = 2
    ),
    "If operator is 'custom', f must be provided"
  )
})
