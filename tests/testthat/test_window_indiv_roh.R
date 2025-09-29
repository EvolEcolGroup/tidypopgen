skip_if_not_installed("detectRUNS")

test_that("windows_indiv_roh works without backingfile update", {
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c"),
    population = c("pop1", "pop1", "pop2")
  )
  test_genotypes <- rbind(
    c(1, 1, 1, 1, 1, 1, 0),
    c(2, 1, 0, 0, 0, 0, 1),
    c(2, 2, 0, 0, 1, 1, 2)
  )
  test_loci <- data.frame(
    name = paste0("rs", 1:7),
    chromosome = paste0("chr", c(1, 1, 1, 1, 2, 2, 2)),
    position = as.integer(c(1000, 1500, 2000, 2500, 50000, 60000, 60500)),
    genetic_dist = as.double(rep(0, 7)),
    allele_ref = c("A", "T", "C", "G", "C", "T", "T"),
    allele_alt = c("T", "C", NA, "C", "G", "A", "A")
  )

  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  )

  suppressMessages(test_roh1 <- windows_indiv_roh(
    test_gt,
    window_size = 4,
    min_snp = 2,
    min_density = 1 / 500,
    max_gap = 5000,
    min_length_bps = 1000,
    threshold = 0.05
  ))

  # test_gt is not grouped, so id and group are identical
  expect_true("group" %in% colnames(test_roh1))
  expect_equal(test_roh1$group, test_roh1$id)

  # remove one snp
  test_gt <- test_gt %>% select_loci(-7)

  # Check roh still works
  suppressMessages(test_roh2 <- windows_indiv_roh(
    test_gt,
    window_size = 4,
    min_snp = 2,
    min_density = 1 / 500,
    max_gap = 5000,
    min_length_bps = 1000,
    threshold = 0.05
  ))

  expect_equal(test_roh1, test_roh2)

  # check that we get deprecation error if we use gt_roh_window
  suppressMessages(expect_warning(
    gt_roh_window(
      test_gt,
      window_size = 4,
      min_snp = 2,
      min_density = 1 / 500,
      max_gap = 5000,
      min_length_bps = 1000,
      threshold = 0.05
    ),
    "This is a soft-deprecated function,"
  ))
})

test_that("grouping method", {
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c"),
    population = c("pop1", "pop1", "pop2")
  )
  test_genotypes <- rbind(
    c(1, 1, 1, 1, 1, 1, 0),
    c(2, 1, 0, 0, 0, 0, 1),
    c(2, 2, 0, 0, 1, 1, 2)
  )
  test_loci <- data.frame(
    name = paste0("rs", 1:7),
    chromosome = paste0("chr", c(1, 1, 1, 1, 2, 2, 2)),
    position = as.integer(c(1000, 1500, 2000, 2500, 50000, 60000, 60500)),
    genetic_dist = as.double(rep(0, 7)),
    allele_ref = c("A", "T", "C", "G", "C", "T", "T"),
    allele_alt = c("T", "C", NA, "C", "G", "A", "A")
  )

  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  )

  test_gt <- test_gt %>% group_by(population)

  suppressMessages(test_roh1 <- windows_indiv_roh(
    test_gt,
    window_size = 4,
    min_snp = 2,
    min_density = 1 / 500,
    max_gap = 5000,
    min_length_bps = 1000,
    threshold = 0.05
  ))

  # test_gt is grouped, so group is grouping levels
  expect_true("group" %in% colnames(test_roh1))
  expect_equal(test_roh1$group, dplyr::group_keys(test_gt) %>% pull(1))
})

test_that("results after subsetting", {
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c", "d"),
    population = c("pop1", "pop1", "pop2", "pop2")
  )
  test_genotypes <- rbind(
    c(1, 1, 1, 1, 1, 1, 0),
    c(2, 1, 0, 0, 0, 0, 1),
    c(2, 2, 2, 1, 0, 1, 2),
    c(2, 2, 0, 0, 1, 1, 2)
  )
  test_loci <- data.frame(
    name = paste0("rs", 1:7),
    chromosome = paste0("chr", c(1, 1, 1, 1, 2, 2, 2)),
    position = as.integer(c(1000, 1500, 2000, 2500, 50000, 60000, 60500)),
    genetic_dist = as.double(rep(0, 7)),
    allele_ref = c("A", "T", "C", "G", "C", "T", "T"),
    allele_alt = c("T", "C", NA, "C", "G", "A", "A")
  )

  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  )

  suppressMessages(test_roh_full <- windows_indiv_roh(
    test_gt,
    window_size = 4,
    min_snp = 2,
    min_density = 1 / 500,
    max_gap = 5000,
    min_length_bps = 1000,
    threshold = 0.05
  ))

  # subset individuals
  test_gt_sub_ind <- test_gt %>% filter(id %in% c("c", "d"))

  suppressMessages(test_roh_sub_ind <- windows_indiv_roh(
    test_gt_sub_ind,
    window_size = 4,
    min_snp = 2,
    min_density = 1 / 500,
    max_gap = 5000,
    min_length_bps = 1000,
    threshold = 0.05
  ))

  # results should be identical to the corresponding individuals in full results
  expect_true(all.equal(test_roh_full[c(2, 3), ],
    test_roh_sub_ind,
    check.attributes = FALSE
  ))
})
