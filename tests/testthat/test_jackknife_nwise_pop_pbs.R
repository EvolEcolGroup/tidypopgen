test_that("jackknife_nwise_pop_pbs matches pbs from genome-wide fst", {
  test_gt <- make_jack_test_gt()
  jack_pbs <- jackknife_nwise_pop_pbs(test_gt,
    type = "tidy", return_fst = TRUE, block_size = 2, size_unit = "snp"
  )
  # pbs from the genome-wide fst
  fst_all <- pairwise_pop_fst(test_gt, method = "Hudson")
  fst_all <- matrix(fst_all$value,
    nrow = 1,
    dimnames = list(NULL, paste0(
      "fst_", fst_all$population_1, ".", fst_all$population_2
    ))
  )
  pbs_all <- unlist(lapply(
    pbs_one_triplet(c("pop1", "pop2", "pop3"), fst_all),
    unname
  ))
  expect_equal(
    jack_pbs$estimate,
    unname(c(fst_all[1, ], pbs_all))
  )
  expect_equal(jack_pbs$stat_name, c(colnames(fst_all), names(pbs_all)))

  # leave-one-out pbs is the same as pbs computed on the loci outside the block
  block_id <- .jackknife_block_ids(show_loci(test_gt),
    block_size = 2, size_unit = "snp"
  )
  jack_matrix <- jackknife_nwise_pop_pbs(test_gt,
    block_size = 2, size_unit = "snp"
  )
  fst_sub <- test_gt %>%
    select_loci(which(block_id != 1)) %>%
    pairwise_pop_fst(method = "Hudson")
  fst_sub <- matrix(fst_sub$value, nrow = 1, dimnames = dimnames(fst_all))
  pbs_sub <- unlist(lapply(
    pbs_one_triplet(c("pop1", "pop2", "pop3"), fst_sub),
    unname
  ))
  fst_blocks <- .jackknife_fst_blocks(test_gt, "Hudson", block_id)
  pbs_loo <- unlist(lapply(
    pbs_one_triplet(
      c("pop1", "pop2", "pop3"),
      fst_blocks$loo[1, , drop = FALSE]
    ),
    unname
  ))
  expect_equal(pbs_loo, pbs_sub)

  # matrix format
  expect_equal(names(jack_matrix), c("quantity", names(pbs_all)))
  expect_equal(unlist(jack_matrix[1, -1]), pbs_all)
})

test_that("jackknife_nwise_pop_pbs checks its inputs", {
  test_gt <- make_jack_test_gt()
  expect_error(
    test_gt %>%
      dplyr::ungroup() %>%
      jackknife_nwise_pop_pbs(block_size = 2),
    ".x should be a grouped gen_tibble"
  )
  test_gt_2pop <- test_gt %>% dplyr::filter(population %in% c("pop1", "pop2"))
  expect_error(
    jackknife_nwise_pop_pbs(test_gt_2pop, block_size = 2),
    "At least 3 populations are required to compute PBS."
  )
})
