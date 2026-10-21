test_that("windows_stats_generic is correct", {
  x <- c(1, 2, 3, 4, 5, 6, 10, 11, 12, 13, 14, 15, 16)
  loci_table <- data.frame(
    name = paste0("locus", 1:13),
    chromosome = c(rep("chr1", 6), rep("chr2", 7)),
    genetic_dist = rep(0, 13),
    position = c(
      50, 120, 150, 180, 230, 390, 110, 120, 150,
      180, 230, 280, 350
    ),
    allele_ref = rep("A", 13),
    allele_alt = rep("T", 13)
  )
  window_res <- windows_stats_generic(
    .x = x,
    loci_table = loci_table,
    window_size = 4,
    step_size = 3,
    size_unit = "snp",
    operator = "sum",
    min_loci = 1
  )
  expect_true(all(window_res$n_loci == c(4, 3, 4, 4)))
  expect_true(window_res$stat[1] == sum(x[1:4]))
  expect_true(window_res$stat[4] == sum(x[10:13]))

  # now set to NA incomplete windows
  window_res <- windows_stats_generic(
    .x = x,
    loci_table = loci_table,
    window_size = 4,
    step_size = 3,
    size_unit = "snp",
    operator = "sum",
    min_loci = 1,
    complete = TRUE
  )
  expect_true((is.na(window_res$stat[2])))

  # now set units to bp
  window_res <- windows_stats_generic(
    .x = x,
    loci_table = loci_table,
    window_size = 100,
    step_size = 50,
    size_unit = "bp",
    operator = "sum",
    min_loci = 1
  )
  # window for chr1 between 251 and 350 is empty (i.e NA)
  expect_true(
    all(
      is.na(
        window_res$stat[
          window_res$chromosome == "chr1" &
            window_res$start == 250
        ]
      )
    )
  )
  # window in chr2 starting at 100 should have n_loci = 4
  expect_true(
    all(
      window_res$n_loci[
        window_res$chromosome == "chr2" &
          window_res$start == 100
      ] == 4
    )
  )
  # the smallest start for chr2 is 100
  expect_true(
    min(
      window_res$start[window_res$chromosome == "chr2"]
    ) == 101
  )
  # window for chr2 and start 251 should have sum 31 (15+16)
  expect_true(
    window_res$stat[
      window_res$chromosome == "chr2" &
        window_res$start == 251
    ] == 31
  )
})

test_that("windows_stats_generic windows include the last loci", {
  x <- rep(1, 100)
  loci_table <- data.frame(
    chromosome = "chr1",
    position = 1:100 * 10
  )
  # (100 - 10) is not a multiple of the step size
  window_res <- windows_stats_generic(
    .x = x,
    loci_table = loci_table,
    window_size = 10,
    step_size = 4,
    size_unit = "snp",
    operator = "sum"
  )
  in_window <- sapply(1:100, function(i) {
    any(window_res$start <= i & window_res$end >= i)
  })
  expect_true(all(in_window))
  # the last window is partial, and includes loci 93 to 100
  expect_equal(utils::tail(window_res$end, 1), 102)
  expect_equal(utils::tail(window_res$n_loci, 1), 8)
  expect_equal(utils::tail(window_res$stat, 1), 8)
  # windows before the last one are unchanged
  expect_equal(window_res$end[-nrow(window_res)], seq(10, 98, by = 4))
  # and it is NA if we only want complete windows
  window_res_complete <- windows_stats_generic(
    .x = x,
    loci_table = loci_table,
    window_size = 10,
    step_size = 4,
    size_unit = "snp",
    operator = "sum",
    complete = TRUE
  )
  expect_true(is.na(utils::tail(window_res_complete$stat, 1)))
  # same in bp
  window_res <- windows_stats_generic(
    .x = x,
    loci_table = loci_table,
    window_size = 100,
    step_size = 40,
    size_unit = "bp",
    operator = "sum"
  )
  in_window <- sapply(loci_table$position, function(i) {
    any(window_res$start <= i & window_res$end >= i)
  })
  expect_true(all(in_window))
  expect_equal(utils::tail(window_res$n_loci, 1), 8)
  # no extra window if the trailing loci are in a gap between windows
  window_res <- windows_stats_generic(
    .x = x,
    loci_table = loci_table,
    window_size = 10,
    step_size = 20,
    size_unit = "snp",
    operator = "sum"
  )
  expect_equal(window_res$end, seq(10, 90, by = 20))
})

test_that("windows generic correctly catches errors", {
  x <- c(1, 2, 3, 4, 5, 6, 10, 11, 12, 13, 14, 15, 16)
  loci_table <- data.frame(
    name = paste0("locus", 1:13),
    chromosome = c(rep("chr1", 6), rep("chr2", 7)),
    genetic_dist = rep(0, 13),
    position = c(
      50, 120, 150, 180, 230, 390, 110, 120, 150,
      180, 230, 280, 350
    ),
    allele_ref = rep("A", 13),
    allele_alt = rep("T", 13)
  )

  expect_error(
    windows_stats_generic(
      .x = x,
      loci_table = loci_table,
      window_size = -1,
      step_size = -1,
      size_unit = "snp",
      operator = "sum",
      min_loci = -1
    ),
    "window_size must be positive"
  )

  expect_error(
    windows_stats_generic(
      .x = x,
      loci_table = loci_table,
      window_size = 4,
      step_size = -1,
      size_unit = "snp",
      operator = "sum",
      min_loci = -1
    ),
    "step_size must be positive"
  )

  expect_error(
    windows_stats_generic(
      .x = x,
      loci_table = loci_table,
      window_size = 4,
      step_size = 1,
      size_unit = "snp",
      operator = "sum",
      min_loci = -1
    ),
    "min_loci must be positive"
  )
  expect_error(
    windows_stats_generic(
      .x = x,
      loci_table = loci_table,
      window_size = 4,
      step_size = 1,
      size_unit = "snp",
      operator = "sum",
      min_loci = 1,
      complete = "blah"
    ),
    "complete must be a boolean"
  )

  # do not test for specific error message from match.arg
  expect_error(
    windows_stats_generic(
      .x = x,
      loci_table = loci_table,
      window_size = 4,
      step_size = 1,
      size_unit = "snp",
      min_loci = 1,
      complete = TRUE,
      operator = "not_a_function"
    )
  )

  # now test if we give incorrectly formatted loci table
  expect_error(
    windows_stats_generic(
      .x = x,
      loci_table = loci_table %>% dplyr::select(-chromosome),
      window_size = 4,
      step_size = 1,
      size_unit = "snp",
      operator = "sum",
      min_loci = 1
    ),
    "loci_table must contain column 'chromosome' when size_unit is 'snp'"
  )

  expect_no_error(
    windows_stats_generic(
      .x = x,
      loci_table = loci_table %>% dplyr::select(-position),
      window_size = 4,
      step_size = 1,
      size_unit = "snp",
      operator = "sum",
      min_loci = 1
    )
  )

  expect_error(
    windows_stats_generic(
      .x = x,
      loci_table = loci_table %>% dplyr::select(-position),
      window_size = 4,
      step_size = 1,
      size_unit = "bp",
      operator = "sum",
      min_loci = 1
    ),
    "loci_table must contain columns 'chromosome' and 'position' "
  )

  expect_error(
    windows_stats_generic(
      .x = x,
      loci_table = loci_table %>% dplyr::select(-position),
      window_size = 2,
      step_size = 1,
      size_unit = "snp",
      operator = "sum",
      min_loci = 3
    ),
    "min_loci must be less than window_size."
  )
})
