# Functions applied to a gen_tibble without individuals, as happens to an
# empty group (from group_by(.drop = FALSE), or filter()/slice() with
# .preserve = TRUE) in group_map(). Their results match those that the
# grouped methods give for an empty group.

test_gt <- load_example_gt("gen_tbl")
n_loci <- count_loci(test_gt)
empty_gt <- test_gt %>% filter(population == "no_such_population")

test_that("loci_* functions give NA for each locus without individuals", {
  expect_equal(nrow(empty_gt), 0)
  freq <- loci_alt_freq(empty_gt)
  expect_length(freq, n_loci)
  expect_true(all(is.na(freq)))
  maf <- loci_maf(empty_gt)
  expect_length(maf, n_loci)
  expect_true(all(is.na(maf)))
  counts <- loci_alt_freq(empty_gt, as_counts = TRUE)
  expect_equal(dim(counts), c(n_loci, 2))
  expect_true(all(counts == 0))
  expect_equal(loci_missingness(empty_gt), rep(NaN, n_loci))
  expect_equal(loci_missingness(empty_gt, as_counts = TRUE), rep(0, n_loci))
  expect_equal(loci_hwe(empty_gt), rep(NA_real_, n_loci))
  pi <- loci_pi(empty_gt)
  expect_length(pi, n_loci)
  expect_true(all(is.na(pi)))
})

test_that("pairwise functions give empty results without individuals", {
  pairwise_funs <- list(
    pairwise_king = pairwise_king,
    pairwise_ibs = pairwise_ibs,
    pairwise_allele_sharing = pairwise_allele_sharing
  )
  for (fun in pairwise_funs) {
    res_mat <- fun(empty_gt, as_matrix = TRUE)
    expect_s3_class(res_mat, "pairwise_matrix")
    expect_equal(dim(res_mat), c(0, 0))
    res_tidy <- fun(empty_gt)
    expect_s3_class(res_tidy, "pairwise_tbl")
    expect_named(res_tidy, c("item1", "item2", "value"))
    expect_equal(nrow(res_tidy), 0)
  }
  expect_equal(dim(pairwise_grm(empty_gt)), c(0, 0))
  expect_error(
    pairwise_ibs(empty_gt, type = "raw_counts"),
    "needs at least one individual"
  )
  expect_length(indiv_inbreeding(empty_gt), 0)
})

test_that("group_map() and grouped indiv_inbreeding() work with empty groups", {
  test_gt$population <- factor(
    test_gt$population,
    levels = c("pop1", "pop_empty", "pop2", "pop3")
  )
  keep_gt <- test_gt %>% group_by(population, .drop = FALSE)
  drop_gt <- test_gt %>% group_by(population)
  # indiv_inbreeding() has an empty result for the empty group
  inb_keep <- indiv_inbreeding(keep_gt)
  expect_length(inb_keep, 4)
  expect_length(inb_keep[[2]], 0)
  expect_equal(inb_keep[-2], indiv_inbreeding(drop_gt))
  # group_map() of loci_* functions: the other groups are unchanged, and the
  # empty group is NA, as in the grouped methods
  for (fun in list(loci_alt_freq, loci_maf, loci_missingness, loci_hwe)) {
    map_keep <- group_map(keep_gt, ~ fun(.x))
    expect_equal(map_keep[-2], group_map(drop_gt, ~ fun(.x)))
    expect_true(all(is.na(map_keep[[2]])))
    grouped_empty <- unclass(fun(keep_gt, type = "matrix"))[, "pop_empty"]
    expect_equal(is.na(map_keep[[2]]), is.na(unname(grouped_empty)))
  }
  # and the pairwise functions give an empty matrix
  king_keep <- group_map(keep_gt, ~ pairwise_king(.x, as_matrix = TRUE))
  expect_equal(dim(king_keep[[2]]), c(0, 0))
})
