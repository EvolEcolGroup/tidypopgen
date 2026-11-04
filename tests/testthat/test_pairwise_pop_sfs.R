# Tests for pairwise_pop_sfs(), against the reference implementation of
# easySFS/dadi in helper_sfs_reference.R

sim <- ref_sfs_simulate(c(popA = 6, popB = 5, popC = 4),
  n_loci = 80,
  seed = 11
)
pops <- c("popA", "popB", "popC")

test_that("pairwise_pop_sfs matches easySFS (folded and unfolded)", {
  pairs <- utils::combn(3, 2)
  for (n_proj in list(c(8, 6, 4), c(12, 10, 8), c(3, 2, 5))) {
    for (folded in c(TRUE, FALSE)) {
      sfs <- pairwise_pop_sfs(sim$gt,
        n_proj = n_proj, folded = folded,
        type = "array"
      )
      expect_named(sfs, c("popA.popB", "popA.popC", "popB.popC"))
      for (i in seq_len(ncol(pairs))) {
        pair <- pairs[, i]
        ref <- ref_sfs(sim$alt_copies, sim$ploidy, sim$population,
          pops = pops[pair], proj = n_proj[pair], polarized = !folded
        )
        expect_equal(sfs_values(sfs[[i]]), sfs_values(ref), tolerance = 1e-10)
        expect_identical(attr(sfs[[i]], "mask"), unname(attr(ref, "mask")))
        expect_identical(dim(sfs[[i]]), as.integer(n_proj[pair] + 1))
        expect_identical(names(dimnames(sfs[[i]])), pops[pair])
        expect_equal(sum(sfs[[i]]), attr(sfs[[i]], "n_loci"))
      }
    }
  }
})

test_that("pairwise_pop_sfs returns a tidy tibble", {
  sfs_tidy <- pairwise_pop_sfs(sim$gt, n_proj = c(4, 3, 2))
  expect_named(sfs_tidy, c(
    "population_1", "population_2", "allele_count_1",
    "allele_count_2", "n_sites", "masked"
  ))
  expect_identical(
    sfs_tidy$population_1,
    rep(c("popA", "popA", "popB"), c(5 * 4, 5 * 3, 4 * 3))
  )
  expect_identical(
    sfs_tidy$population_2,
    rep(c("popB", "popC", "popC"), c(5 * 4, 5 * 3, 4 * 3))
  )
  sfs_array <- pairwise_pop_sfs(sim$gt, n_proj = c(4, 3, 2), type = "array")
  # by default, the order of the populations is kept: rows are sorted by
  # the first population, so the second changes fastest
  ab <- sfs_tidy[sfs_tidy$population_2 == "popB", ]
  expect_identical(ab$allele_count_1, rep(0:4, each = 4))
  expect_identical(ab$allele_count_2, rep(0:3, times = 5))
  expect_equal(ab$n_sites, as.vector(t(sfs_array$popA.popB)))
  expect_equal(ab$masked, as.vector(t(attr(sfs_array$popA.popB, "mask"))))
  # or the first population changes fastest
  sfs_first <- pairwise_pop_sfs(sim$gt, n_proj = c(4, 3, 2), fastest = "first")
  ab_first <- sfs_first[sfs_first$population_2 == "popB", ]
  expect_identical(ab_first$allele_count_1, rep(0:4, times = 4))
  expect_identical(ab_first$allele_count_2, rep(0:3, each = 5))
  expect_equal(ab_first$n_sites, as.vector(sfs_array$popA.popB))
  expect_identical(attr(sfs_tidy, "n_proj"), c(popA = 4L, popB = 3L, popC = 2L))
})

test_that("pairwise_pop_sfs marginals match pop_sfs without missing data", {
  sim_full <- ref_sfs_simulate(c(popA = 4, popB = 3),
    n_loci = 40,
    max_missing = 0, seed = 12
  )
  sfs_2d <- pairwise_pop_sfs(sim_full$gt,
    n_proj = c(6, 4), folded = FALSE,
    type = "array"
  )[[1]]
  sfs_1d <- pop_sfs(sim_full$gt,
    n_proj = c(6, 4), folded = FALSE,
    type = "array"
  )
  expect_equal(rowSums(sfs_2d), as.vector(sfs_1d$popA), ignore_attr = TRUE)
  expect_equal(colSums(sfs_2d), as.vector(sfs_1d$popB), ignore_attr = TRUE)
})

test_that("pairwise_pop_sfs does not depend on the number of cores", {
  sfs_1 <- pairwise_pop_sfs(sim$gt,
    n_proj = c(8, 6, 4), n_cores = 1,
    type = "array"
  )
  sfs_2 <- pairwise_pop_sfs(sim$gt,
    n_proj = c(8, 6, 4), n_cores = 2,
    block_size = 9, type = "array"
  )
  expect_identical(sfs_1, sfs_2)
})

test_that("pairwise_pop_sfs chooses the projections by default", {
  best <- pop_sfs_preview(sim$gt) %>% dplyr::filter(best)
  expect_message(
    sfs_default <- pairwise_pop_sfs(sim$gt),
    "maximises its number of segregating sites"
  )
  expect_identical(
    sfs_default,
    pairwise_pop_sfs(sim$gt, n_proj = stats::setNames(best$n_proj, best$group))
  )
})

test_that("pairwise_pop_sfs checks its arguments", {
  expect_error(
    pairwise_pop_sfs(dplyr::ungroup(sim$gt), n_proj = 2),
    "grouped gen_tibble"
  )
  expect_error(
    pairwise_pop_sfs(sim$gt %>% dplyr::filter(population == "popA"),
      n_proj = 2
    ),
    "at least 2 groups"
  )
  expect_error(pairwise_pop_sfs(sim$gt, n_proj = 13), "popA \\(12\\)")
})

test_that("pairwise_pop_sfs uses one column per grouping variable", {
  gt_two <- sim$gt %>%
    dplyr::ungroup() %>%
    dplyr::mutate(sex = rep(c("F", "M"), length.out = dplyr::n())) %>%
    dplyr::group_by(population, sex)
  sfs <- pairwise_pop_sfs(gt_two, n_proj = 2)
  expect_named(sfs, c(
    "population_1", "sex_1", "population_2", "sex_2",
    "allele_count_1", "allele_count_2", "n_sites", "masked"
  ))
  expect_identical(nrow(sfs), as.integer(choose(6, 2) * 9))
  sfs_array <- pairwise_pop_sfs(gt_two, n_proj = 2, type = "array")
  expect_identical(names(sfs_array)[1], "popA.F.popA.M")
})
