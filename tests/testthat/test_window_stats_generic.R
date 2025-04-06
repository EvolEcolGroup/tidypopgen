test_that("window_stats_generic is correct", {
  x <- c(1, 2, 3, 4, 5, 6, 10, 11, 12, 13, 14, 15, 16)
  loci_table <- data.frame(
    name = paste0("locus", 1:13),
    chromosome = c(rep("chr1", 6), rep("chr2", 7)),
    genetic_dist = rep(0, 13),
    position = c(50, 120, 150, 180, 230, 390, 110, 120, 150,
                 180, 230, 280, 350),
    allele_ref = rep("A", 13),
    allele_alt = rep("T", 13)
  )
  window_res <- window_stats_generic(
    x = x,
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
  window_res <- window_stats_generic(
    x = x,
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
  window_res <- window_stats_generic(
    x = x,
    loci_table = loci_table,
    window_size = 100,
    step_size = 50,
    size_unit = "bp",
    operator = "sum",
    min_loci = 1
  )
  # window for chr1 between 251 and 350 is empty (i.e NA)
  expect_true(all(is.na(window_res$stat[window_res$chromosome == "chr1" &
                                          window_res$start == 250])))
  # window in chr2 starting at 100 should have n_loci = 4
  expect_true(all(window_res$n_loci[window_res$chromosome == "chr2" &
                                      window_res$start == 100] == 4))
  # the smallest start for chr2 is 100
  expect_true(min(window_res$start[window_res$chromosome == "chr2"]) == 101)
  # window for chr2 and start 251 should have sum 31 (15+16)
  expect_true(window_res$stat[window_res$chromosome == "chr2" &
                                window_res$start == 251] == 31)
})
