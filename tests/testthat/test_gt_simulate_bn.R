test_that("gt_simulate_bn gives populations with the requested structure", {
  set.seed(123)
  sim_gt <- gt_simulate_bn(
    n_indiv = c(30, 40), n_loci = 2000, fst = c(0.01, 0.2)
  )
  expect_s3_class(sim_gt, "gen_tbl")
  expect_equal(nrow(sim_gt), 70)
  expect_equal(count_loci(sim_gt), 2000)
  expect_equal(as.vector(table(sim_gt$population)), c(30, 40))
  expect_false(anyNA(show_genotypes(sim_gt)))
  # the more drifted population is further from the ancestral frequencies,
  # so it has lower heterozygosity (expected 2p(1-p)(1-F))
  het <- sim_gt %>%
    dplyr::group_by(population) %>%
    pop_het_obs()
  expect_gt(het[1], het[2])
  # the same seed gives the same dataset
  set.seed(123)
  again <- gt_simulate_bn(
    n_indiv = c(30, 40), n_loci = 2000, fst = c(0.01, 0.2)
  )
  expect_equal(show_genotypes(again), show_genotypes(sim_gt))
})

test_that("gt_simulate_bn sets missing genotypes per individual", {
  set.seed(42)
  rate <- rep(c(0, 0.5, 0.9), each = 10)
  sim_gt <- gt_simulate_bn(n_indiv = 30, n_loci = 3000, fst = 0, missing = rate)
  observed <- indiv_missingness(sim_gt)
  expect_equal(observed[1:10], rep(0, 10))
  expect_equal(mean(observed[11:20]), 0.5, tolerance = 0.02)
  expect_equal(mean(observed[21:30]), 0.9, tolerance = 0.02)
})

test_that("gt_simulate_bn checks its arguments", {
  expect_error(
    gt_simulate_bn(n_indiv = c(10, 10), fst = c(0.1, 0.1, 0.1)),
    "one per population"
  )
  expect_error(gt_simulate_bn(fst = 1), "'fst' should be")
  expect_error(
    gt_simulate_bn(n_indiv = 10, missing = c(0.1, 0.2)),
    "one per individual"
  )
  expect_error(gt_simulate_bn(freq_range = c(0.9, 0.1)), "lower first")
})
