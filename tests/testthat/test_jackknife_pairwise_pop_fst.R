test_that("jackknife_pairwise_pop_fst leave-one-out matches subsets", {
  test_gt <- make_jack_test_gt()
  for (method in c("Hudson", "Nei87", "WC84")) {
    block_id <- .jackknife_block_ids(show_loci(test_gt),
      block_size = 1, size_unit = "cM"
    )
    fst_blocks <- .jackknife_fst_blocks(test_gt,
      method = method, block_id = block_id
    )
    # the full estimate is the genome-wide fst
    fst_all <- pairwise_pop_fst(test_gt, method = method)
    expect_equal(unname(fst_blocks$full), fst_all$value)
    # leave-one-out estimates are the same as the genome-wide fst
    # computed on the loci outside the block
    for (i_block in unique(block_id)) {
      fst_sub <- test_gt %>%
        select_loci(which(block_id != i_block)) %>%
        pairwise_pop_fst(method = method)
      expect_equal(unname(fst_blocks$loo[i_block, ]), fst_sub$value)
    }
  }
})

test_that("jackknife_pairwise_pop_fst returns the right format", {
  test_gt <- make_jack_test_gt()
  jack_matrix <- jackknife_pairwise_pop_fst(test_gt,
    method = "WC84", block_size = 2, size_unit = "snp"
  )
  jack_tidy <- jackknife_pairwise_pop_fst(test_gt,
    type = "tidy", method = "WC84", block_size = 2, size_unit = "snp"
  )
  expect_equal(
    names(jack_matrix),
    c("quantity", "fst_pop1.pop2", "fst_pop1.pop3", "fst_pop2.pop3")
  )
  expect_equal(
    jack_matrix$quantity,
    c("estimate", "jack_estimate", "se", "z", "n_blocks")
  )
  expect_equal(
    names(jack_tidy),
    c("stat_name", "estimate", "jack_estimate", "se", "z", "n_blocks")
  )
  expect_equal(jack_tidy$stat_name, names(jack_matrix)[-1])
  expect_equal(jack_tidy$se, unname(unlist(jack_matrix[3, -1])))
  # estimate is the genome-wide fst
  fst_all <- pairwise_pop_fst(test_gt, method = "WC84")
  expect_equal(jack_tidy$estimate, fst_all$value)
  expect_true(all(jack_tidy$n_blocks == 5))
  expect_true(all(jack_tidy$se > 0))
})
