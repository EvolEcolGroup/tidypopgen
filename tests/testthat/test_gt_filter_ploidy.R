# tests for .gt_refresh_ploidy(), used by filter.gen_tbl() and
# filter.grouped_gen_tbl() to keep the top-level ploidy label from going
# stale once individuals are dropped from a mixed-ploidy or pseudohaploid
# gen_tibble.

mixed_ploidy_gt <- function() {
  # a,b diploid; c,d tetraploid
  test_genotypes <- rbind(
    c(1, 2, 0, 1),
    c(0, 2, 1, 0),
    c(4, 4, 0, 2),
    c(0, 1, 4, 4)
  )
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c", "d"),
    population = c("pop1", "pop1", "pop2", "pop2")
  )
  test_loci <- data.frame(
    name = paste0("rs", 1:4),
    chromosome = c(1, 1, 2, 2),
    position = c(3, 5, 23, 456),
    genetic_dist = as.double(rep(0, 4)),
    allele_ref = c("A", "T", "C", "G"),
    allele_alt = c("T", "C", "G", "A")
  )
  gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    ploidy = c(2, 2, 4, 4),
    quiet = TRUE
  )
}

test_that("filter refreshes ploidy when a mixed-ploidy gen_tibble becomes uniform", { # nolint
  test_gt <- mixed_ploidy_gt()
  expect_equal(show_ploidy(test_gt), 0)

  # keep only the diploid individuals
  dip_gt <- test_gt %>% filter(id %in% c("a", "b"))
  expect_equal(show_ploidy(dip_gt), 2)
  expect_equal(indiv_ploidy(dip_gt), c(2, 2))
  expect_true(tidypopgen:::is_diploid_only(dip_gt))

  # a function gated by stopifnot_diploid() now works, purely because the
  # label was refreshed to a definite 2 (it would otherwise still be
  # blocked by the stale 0, even though the data are now all diploid)
  expect_no_error(indiv_inbreeding(dip_gt))

  # keep only the tetraploid individuals: refreshed to 4, not hardcoded to 2
  tetra_gt <- test_gt %>% filter(id %in% c("c", "d"))
  expect_equal(show_ploidy(tetra_gt), 4)
  expect_equal(indiv_ploidy(tetra_gt), c(4, 4))
})

test_that("filter leaves ploidy at 0 when the remaining individuals are still mixed", { # nolint
  test_gt <- mixed_ploidy_gt()

  # drop one diploid and one tetraploid individual: still mixed
  still_mixed_gt <- test_gt %>% filter(id %in% c("a", "c"))
  expect_equal(show_ploidy(still_mixed_gt), 0)
  expect_equal(indiv_ploidy(still_mixed_gt), c(2, 4))
})

test_that("filter refreshes ploidy for grouped mixed-ploidy gen_tibbles", {
  test_gt <- mixed_ploidy_gt() %>% dplyr::group_by(population)
  expect_equal(show_ploidy(test_gt), 0)

  # pop1 is entirely diploid (a, b); pop2 is entirely tetraploid (c, d)
  dip_gt <- test_gt %>% filter(population == "pop1")
  expect_equal(show_ploidy(dip_gt), 2)

  tetra_gt <- test_gt %>% filter(population == "pop2")
  expect_equal(show_ploidy(tetra_gt), 4)
})

test_that("filter refreshes ploidy for grouped pseudohaploid gen_tibbles", {
  # filter.grouped_gen_tbl() previously did not refresh the ploidy label at
  # all (only filter.gen_tbl() did), so a grouped filter that dropped all
  # pseudohaploid individuals left the gen_tibble incorrectly flagged -2
  test_genotypes <- rbind(
    c(2, NA, 0, NA, 2, 0),
    c(2, 2, 0, NA, 0, NA),
    c(2, NA, 0, 0, 1, 1),
    c(2, 0, 0, 1, 0, 0),
    c(1, 2, 0, 1, 2, 1),
    c(0, 0, 1, 0, 1, 2),
    c(0, 1, 2, 0, 1, 1)
  )
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c", "d", "e", "f", "g"),
    population = c("pop1", "pop1", "pop2", "pop2", "pop1", "pop3", "pop3")
  )
  test_loci <- data.frame(
    name = paste0("rs", 1:6),
    chromosome = paste0("chr", c(1, 1, 1, 1, 2, 2)),
    position = as.integer(c(3, 5, 65, 343, 23, 456)),
    genetic_dist = as.double(rep(0, 6)),
    allele_ref = c("A", "T", "C", "G", "C", "T"),
    allele_alt = c("T", "C", NA, "C", "G", "A")
  )
  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  )
  test_gt <- gt_pseudohaploid(test_gt) %>% dplyr::group_by(population)
  expect_equal(show_ploidy(test_gt), -2)

  clean_gt <- test_gt %>% filter(indiv_missingness(genotypes) < 0.15)
  expect_equal(show_ploidy(clean_gt), 2)
  expect_no_error(indiv_het_obs(clean_gt))
})

test_that("filter does not error when it drops all individuals", {
  test_gt <- mixed_ploidy_gt()
  empty_gt <- test_gt %>% filter(id == "nonexistent")
  expect_equal(nrow(empty_gt), 0)
  # ploidy label is left untouched (still mixed) rather than erroring
  expect_equal(show_ploidy(empty_gt), 0)
})
