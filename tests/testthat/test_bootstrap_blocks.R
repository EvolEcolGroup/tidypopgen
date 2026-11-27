test_that("bootstrap_blocks draws blocks that never span chromosomes", {
  loci <- data.frame(
    chromosome = c(1, 1, 1, 2, 2, 2, 2),
    position = c(10, 20, 150, 5, 60, 120, 250)
  )
  set.seed(1)
  blocks <- bootstrap_blocks(loci,
    n_bootstrap = 4, block_size = 2,
    size_unit = "snp"
  )
  # blocks of 2 loci within each chromosome, numbered in order
  expect_identical(blocks$block_id, c(1L, 1L, 2L, 3L, 3L, 4L, 4L))
  # each replicate draws as many blocks as there are
  expect_identical(dim(blocks$weights), c(4L, 4L))
  expect_true(is.integer(blocks$weights))
  expect_identical(rowSums(blocks$weights), rep(4, 4))
  # the same blocks after the same seed
  set.seed(1)
  expect_identical(
    bootstrap_blocks(loci, n_bootstrap = 4, block_size = 2, size_unit = "snp"),
    blocks
  )
  # blocks in base pairs: positions (0, 100], (100, 200], ... of each
  # chromosome
  bp <- bootstrap_blocks(loci,
    n_bootstrap = 2, block_size = 100,
    size_unit = "bp"
  )
  expect_identical(bp$block_id, c(1L, 1L, 2L, 3L, 3L, 4L, 5L))
  expect_identical(ncol(bp$weights), 5L)
})

test_that("bootstrap_blocks needs at least 2 blocks", {
  loci <- data.frame(chromosome = c(1, 1, 1), position = c(10, 20, 30))
  expect_error(
    bootstrap_blocks(loci, n_bootstrap = 2, block_size = 10, size_unit = "snp"),
    "at least 2 blocks"
  )
  expect_error(
    bootstrap_blocks(loci, n_bootstrap = 2, size_unit = "snp"),
    "`block_size` is needed"
  )
})
