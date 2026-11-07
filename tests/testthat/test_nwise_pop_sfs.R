# Tests for nwise_pop_sfs(), against the reference implementation of
# easySFS/dadi in helper_sfs_reference.R

sim <- ref_sfs_simulate(c(popA = 5, popB = 4, popC = 3),
  n_loci = 80,
  seed = 21
)
pops <- c("popA", "popB", "popC")

test_that("nwise_pop_sfs matches easySFS (folded and unfolded)", {
  for (n_proj in list(c(6, 4, 4), c(10, 8, 6), c(2, 3, 1))) {
    for (folded in c(TRUE, FALSE)) {
      sfs <- nwise_pop_sfs(sim$gt,
        n_proj = n_proj, folded = folded,
        type = "array"
      )
      ref <- ref_sfs(sim$alt_copies, sim$ploidy, sim$population,
        pops = pops, proj = n_proj, polarized = !folded
      )
      expect_s3_class(sfs, "sfs")
      expect_equal(sfs_values(sfs), sfs_values(ref), tolerance = 1e-10)
      expect_identical(attr(sfs, "mask"), unname(attr(ref, "mask")))
      expect_identical(names(dimnames(sfs)), pops)
      expect_equal(sum(sfs), attr(sfs, "n_loci"))
    }
  }
})

test_that("nwise_pop_sfs keeps the order of the populations", {
  sfs_tidy <- nwise_pop_sfs(sim$gt, n_proj = c(3, 2, 1))
  expect_named(sfs_tidy, c(pops, "n_sites", "masked"))
  # rows are sorted by the first population, then the second, ... (the last
  # population changes fastest)
  expect_identical(sfs_tidy$popA, rep(0:3, each = 3 * 2))
  expect_identical(sfs_tidy$popB, rep(rep(0:2, each = 2), times = 4))
  expect_identical(sfs_tidy$popC, rep(0:1, times = 4 * 3))
  sfs_array <- nwise_pop_sfs(sim$gt, n_proj = c(3, 2, 1), type = "array")
  expect_equal(sfs_tidy$n_sites, as.vector(aperm(sfs_array, 3:1)))
  expect_identical(
    sfs_tidy$masked,
    as.vector(aperm(attr(sfs_array, "mask"), 3:1))
  )
  # tidy() of the array is the same
  expect_identical(tidy(sfs_array), tibble::as_tibble(sfs_tidy))
  # the first population can change fastest instead
  sfs_first <- nwise_pop_sfs(sim$gt, n_proj = c(3, 2, 1), fastest = "first")
  expect_identical(sfs_first$popA, rep(0:3, times = 3 * 2))
  expect_identical(sfs_first$popC, rep(0:1, each = 4 * 3))
  expect_equal(sfs_first$n_sites, as.vector(sfs_array))
  expect_identical(tidy(sfs_array, fastest = "first"), tibble::as_tibble(
    sfs_first
  ))
  # reordering the groups permutes the dimensions of the spectrum
  gt_reordered <- sim$gt %>%
    dplyr::ungroup() %>%
    dplyr::mutate(population = factor(population,
      levels = c("popC", "popA", "popB")
    )) %>%
    dplyr::group_by(population)
  sfs_reordered <- nwise_pop_sfs(gt_reordered,
    n_proj = c(popA = 3, popB = 2, popC = 1), type = "array"
  )
  expect_identical(names(dimnames(sfs_reordered)), c("popC", "popA", "popB"))
  expect_equal(
    sfs_values(sfs_reordered),
    sfs_values(aperm(unclass(sfs_array), c(3, 1, 2)))
  )
})

test_that("nwise_pop_sfs of two groups is their pairwise spectrum", {
  gt_two <- sim$gt %>% dplyr::filter(population != "popC")
  expect_identical(
    nwise_pop_sfs(gt_two, n_proj = c(6, 4), type = "array"),
    pairwise_pop_sfs(gt_two, n_proj = c(6, 4), type = "array")[[1]]
  )
})

test_that("nwise_pop_sfs checks its arguments", {
  expect_error(
    nwise_pop_sfs(sim$gt, n_proj = c(10, 8, 6), max_cells = 100),
    "693 cells"
  )
  expect_error(
    nwise_pop_sfs(dplyr::ungroup(sim$gt), n_proj = 2),
    "grouped gen_tibble"
  )
  expect_error(
    nwise_pop_sfs(sim$gt %>% dplyr::filter(population == "popA"), n_proj = 2),
    "at least 2 groups"
  )
  expect_message(
    nwise_pop_sfs(sim$gt),
    "maximises its number of segregating sites"
  )
})

test_that("nwise_pop_sfs works with many populations", {
  # with 6 populations of 50 diploids, the configurations of a locus do not
  # fit in a 64 bit key, so the kernel uses a vector key
  sim_many <- ref_sfs_simulate(
    stats::setNames(rep(50, 6), paste0("pop", 1:6)),
    n_loci = 15, seed = 22
  )
  n_proj <- rep(2, 6)
  sfs <- nwise_pop_sfs(sim_many$gt, n_proj = n_proj, type = "array")
  ref <- ref_sfs(sim_many$alt_copies, sim_many$ploidy, sim_many$population,
    pops = paste0("pop", 1:6), proj = n_proj
  )
  expect_equal(sfs_values(sfs), sfs_values(ref), tolerance = 1e-10)
})

test_that("pop_sfs works with large populations", {
  # with more than about 2000 allele copies, the 1D configurations are
  # counted in a hash table rather than a dense vector
  sim_large <- ref_sfs_simulate(c(popA = 1100), n_loci = 10, seed = 23)
  sfs <- pop_sfs(sim_large$gt, n_proj = 1500, type = "array")
  ref <- ref_sfs(sim_large$alt_copies, sim_large$ploidy,
    sim_large$population,
    pops = "popA", proj = 1500
  )
  expect_equal(sfs_values(sfs[[1]]), sfs_values(ref), tolerance = 1e-8)
})
