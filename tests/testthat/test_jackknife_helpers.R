test_that(".jackknife_block_ids never spans chromosomes", {
  test_gt <- make_jack_test_gt()
  loci <- show_loci(test_gt)
  # snps
  block_snp <- .jackknife_block_ids(loci, block_size = 2, size_unit = "snp")
  expect_equal(block_snp, c(1, 1, 2, 3, 3, 4, 4, 5))
  # bp
  block_bp <- .jackknife_block_ids(loci, block_size = 100, size_unit = "bp")
  expect_equal(block_bp, c(1, 1, 2, 3, 3, 4, 5, 6))
  # cM
  block_cm <- .jackknife_block_ids(loci, block_size = 1, size_unit = "cM")
  expect_equal(block_cm, c(1, 1, 2, 3, 3, 4, 5, 6))
  # no block includes loci from more than one chromosome
  for (block_id in list(block_snp, block_bp, block_cm)) {
    expect_true(all(tapply(loci$chromosome, block_id, function(x) {
      length(unique(x))
    }) == 1))
  }
})

test_that(".jackknife_block_ids checks its inputs", {
  test_gt <- make_jack_test_gt()
  loci <- show_loci(test_gt)
  expect_error(
    .jackknife_block_ids(loci, block_size = -1, size_unit = "snp"),
    "block_size must be a single positive number"
  )
  loci$genetic_dist <- 0
  expect_error(
    .jackknife_block_ids(loci, block_size = 1, size_unit = "cM"),
    "all values of 'genetic_dist' are 0"
  )
  expect_error(
    .jackknife_block_ids(loci[, c("chromosome", "position")],
      block_size = 1, size_unit = "cM"
    ),
    "must contain columns 'chromosome' and 'genetic_dist'"
  )
})

test_that(".jackknife_weighted matches the classic jackknife for equal m", {
  theta_full <- 0.5
  theta_loo <- matrix(c(0.45, 0.52, 0.49, 0.55, 0.47), ncol = 1)
  m <- matrix(10, nrow = 5, ncol = 1)
  g <- 5
  jack <- .jackknife_weighted(theta_full, theta_loo, m)
  # classic delete-1-block jackknife
  classic_est <- g * theta_full - (g - 1) * mean(theta_loo)
  classic_se <- sqrt((g - 1) / g * sum((theta_loo - mean(theta_loo))^2))
  expect_equal(jack$estimate, theta_full)
  expect_equal(jack$jack_estimate, classic_est)
  expect_equal(jack$se, classic_se)
  expect_equal(jack$z, classic_est / classic_se)
  expect_equal(jack$n_blocks, 5L)
})

test_that(".jackknife_weighted follows Busing et al. for unequal m", {
  theta_full <- c(1.2, 3)
  theta_loo <- cbind(c(1.1, 1.3, 1.25), c(2.9, 3.2, 3.1))
  m <- cbind(c(2, 5, 3), c(4, 0, 6))
  jack <- .jackknife_weighted(theta_full, theta_loo, m)
  # first statistic, by hand
  n <- 10
  h <- n / c(2, 5, 3)
  est <- 3 * 1.2 - sum((1 - c(2, 5, 3) / n) * c(1.1, 1.3, 1.25))
  tau <- h * 1.2 - (h - 1) * c(1.1, 1.3, 1.25)
  se <- sqrt(sum((tau - est)^2 / (h - 1)) / 3)
  expect_equal(jack$jack_estimate[1], est)
  expect_equal(jack$se[1], se)
  # second statistic: the empty block is ignored
  n <- 10
  h <- n / c(4, 6)
  est <- 2 * 3 - sum((1 - c(4, 6) / n) * c(2.9, 3.1))
  tau <- h * 3 - (h - 1) * c(2.9, 3.1)
  se <- sqrt(sum((tau - est)^2 / (h - 1)) / 2)
  expect_equal(jack$jack_estimate[2], est)
  expect_equal(jack$se[2], se)
  expect_equal(jack$n_blocks, c(3L, 2L))
  # with a single block, we can not compute a standard error
  jack_one <- .jackknife_weighted(1, matrix(1), matrix(3))
  expect_true(is.na(jack_one$se))
})

test_that(".block_sums does not create row names", {
  test_gt <- make_jack_test_gt()
  block_id <- .jackknife_block_ids(show_loci(test_gt),
    block_size = 2, size_unit = "snp"
  )
  # a named vector of per locus values (e.g. with loci names)
  x <- stats::setNames(
    c(0.1, 0.4, NA, 0.3, 0.2, 0.6, 0.5, 0.1),
    loci_names(test_gt)
  )
  x_blocks <- .block_sums(x, block_id)
  expect_null(rownames(x_blocks$sums))
  expect_null(colnames(x_blocks$sums))
  expect_null(rownames(x_blocks$n))
  expect_equal(as.vector(x_blocks$sums), c(0.5, 0, 0.5, 1.1, 0.1))
  expect_equal(as.vector(x_blocks$n), c(2, 0, 2, 2, 1))
  # only the column names (one per statistic) are kept for matrices
  fst_blocks <- .jackknife_fst_blocks(test_gt, "Hudson", block_id)
  expect_null(rownames(fst_blocks$loo))
  expect_null(rownames(fst_blocks$m))
})
