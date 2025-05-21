test_genotypes <- rbind(
  c(1, 0, 1, 1, 0, 2, 1, 0),
  c(1, 0, NA, 0, 0, 0, 1, 2),
  c(NA, 0, 0, 1, 1, 0, 1, 0),
  c(0, 0, 1, 0, 0, 0, 1, 0),
  c(2, 0, 1, 2, 1, 1, 2, 2),
  c(0, 0, 0, NA, 1, 0, 0, 1),
  c(1, 1, 0, 1, NA, 0, 1, 0),
  c(0, 0, 1, 0, 0, 0, 1, 2),
  c(1, 0, 1, 0, 0, 0, 1, 0),
  c(0, 0, 1, 0, 0, 0, 1, 2),
  c(1, 1, 1, 0, 1, 0, 1, 2),
  c(0, 1, 0, 1, 1, 0, 2, 2)
)
test_loci <- data.frame(
  name = paste0("rs", 1:8),
  chromosome = paste0("chr", c(1, 1, 1, 2, 2, 2, 2, 2)),
  position = as.integer(c(5, 65, 343, 23, 56, 138, 230, 456)),
  genetic_dist = as.double(rep(0, 8)),
  allele_ref = c("T", "C", "G", "C", "T", "A", "C", "G"),
  allele_alt = c("C", NA, "C", "G", "A", "T", "A", "C")
)
test_indiv_meta <- data.frame(
  id = c("a", "b", "c", "d", "e", "f", "g", "h", "i", "l", "m", "n"),
  population = c(
    "pop1", "pop1", "pop2", "pop2", "pop1", "pop3", "pop3",
    "pop2", "pop3", "pop3", "pop4", "pop4"
  )
)
test_gt <- gen_tibble(
  x = test_genotypes,
  loci = test_loci,
  indiv_meta = test_indiv_meta,
  quiet = TRUE
)
test_gt <- test_gt %>% dplyr::group_by(population)

test_that("windows_nwise_pop_pbs works correctly", {
  test_window <- test_gt %>%
    windows_nwise_pop_pbs(
      window_size = 3,
      step_size = 2,
      size_unit = "snp",
      min_loci = 2
    )
  # expect that we have 4 (combinations) * 6 (stats) columns +3 (window info)
  expect_equal(ncol(test_window), 4 * 6 + 3)

  # TODO check that the numbers are correct

  # if we ask to return fst, we expect more columns (6 pairwise fst extra)
  test_window_fst <- test_gt %>%
    windows_nwise_pop_pbs(
      window_size = 3,
      step_size = 2,
      size_unit = "snp",
      min_loci = 2,
      return_fst = TRUE
    )
  expect_true(ncol(test_window_fst) == ncol(test_window) + 6)

  # expect error when using an ungrouped tibble
  test_gt_ungrouped <- test_gt %>% ungroup()
  expect_error(
    test_gt_ungrouped %>%
      windows_nwise_pop_pbs(
        window_size = 3,
        step_size = 2,
        size_unit = "snp",
        min_loci = 2
      ),
    ".x should be a grouped gen_tibble"
  )

  # expect error if we only have two populations
  test_gt_2pop <- test_gt %>% dplyr::filter(population %in% c("pop1", "pop2"))
  expect_error(
    test_gt_2pop %>%
      windows_nwise_pop_pbs(
        window_size = 3,
        step_size = 2,
        size_unit = "snp",
        min_loci = 2
      ),
    "At least 3 populations are required to compute PBS."
  )
})

test_that("windows type", {
  pbs_windows_matrix <- test_gt %>%
    windows_nwise_pop_pbs(
      type = "matrix",
      window_size = 3,
      step_size = 2,
      size_unit = "snp",
      min_loci = 2
    )
  expect_true(is.data.frame(pbs_windows_matrix))
  pbs_windows_tidy <- test_gt %>%
    windows_nwise_pop_pbs(
      type = "tidy",
      window_size = 3,
      step_size = 2,
      size_unit = "snp",
      min_loci = 2
    )
  expect_true(is.data.frame(pbs_windows_tidy))
  # Compare
  pop1_pop2_tidy <-
    subset(pbs_windows_tidy, pbs_windows_tidy$stat_name == "pbs_pop1.pop2.pop3")
  expect_equal(pop1_pop2_tidy$value, pbs_windows_matrix$pbs_pop1.pop2.pop3)

  # Now with return_fst
  pbs_windows_matrix_fst <- test_gt %>%
    windows_nwise_pop_pbs(
      type = "matrix",
      window_size = 3,
      step_size = 2,
      size_unit = "snp",
      min_loci = 2,
      return_fst = TRUE
    )
  expect_true(is.data.frame(pbs_windows_matrix))
  pbs_windows_tidy_fst <- test_gt %>%
    windows_nwise_pop_pbs(
      type = "tidy",
      window_size = 3,
      step_size = 2,
      size_unit = "snp",
      min_loci = 2,
      return_fst = TRUE
    )
  expect_true(is.data.frame(pbs_windows_tidy_fst))
})
