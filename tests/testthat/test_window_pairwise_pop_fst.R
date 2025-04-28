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
  c(0, 0, 1, 0, 0, 0, 1, 2)
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
  id = c("a", "b", "c", "d", "e", "f", "g", "h", "i", "l"),
  population = c(
    "pop1", "pop1", "pop2", "pop2", "pop1", "pop3", "pop3",
    "pop2", "pop3", "pop3"
  )
)
test_gt <- gen_tibble(
  x = test_genotypes,
  loci = test_loci,
  indiv_meta = test_indiv_meta,
  quiet = TRUE
)
test_gt <- test_gt %>% dplyr::group_by(population)


# testing the infrastructure for windowing
test_that("pairwise_pop_fst num and dem are returned correctly", {
  hudson_gt_fst <- test_gt %>%
    pairwise_pop_fst(method = "Hudson", by_locus = TRUE)
  hudson_gt_fst_nd <- test_gt %>%
    pairwise_pop_fst(
      method = "Hudson",
      return_num_dem = TRUE
    )
  # check that the nums and denominators are indeed the right ones to
  # generate the Fst for each locus
  hudson_gt_from_num_dem <- hudson_gt_fst_nd$Fst_by_locus_num /
    hudson_gt_fst_nd$Fst_by_locus_den
  colnames(hudson_gt_from_num_dem) <- paste0("fst_",
                                             colnames(hudson_gt_from_num_dem))
  expect_equal(hudson_gt_from_num_dem, hudson_gt_fst$Fst_by_locus)
  # confirm error if trying to get num and dem with other methods
  expect_error(
    test_gt %>% pairwise_pop_fst(method = "Nei87", return_num_dem = TRUE),
    "return_num_dem is only available for Hudson"
  )
})


test_that("window_pairwise_pop_fst works correctly", {
  test_window <- test_gt %>%
    window_pairwise_pop_fst(
      window_size = 3,
      step_size = 2,
      size_unit = "snp",
      min_loci = 2
    )

  hudson_gt_fst_nd <- test_gt %>%
    pairwise_pop_fst(
      method = "Hudson",
      return_num_dem = TRUE
    )
  # TODO check that this works correctly!

  bp_window <- window_pairwise_pop_fst(test_gt,
    window_size = 200,
    step_size = 100,
    size_unit = "bp",
    min_loci = 1
  )

  # TODO check that this works correctly
})
