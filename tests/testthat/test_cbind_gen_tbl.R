test_that("cbind correctly binds", {
  test_genotypes <- rbind(
    c(1, 1, 0, 1, 1, 0),
    c(2, 1, 0, NA, 0, 0),
    c(2, NA, 0, 0, 1, 1),
    c(1, 0, 0, 1, 0, 0),
    c(1, 2, 0, 1, 2, 1),
    c(0, 0, 0, 0, NA, 1),
    c(0, 1, 1, 0, 1, NA)
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
  # new data.frame to cbind
  new_df <- data.frame(
    id2 = c("a", "b", "c", "d", "e", "f", "g"),
    population = c(
      "region1", "region1", "region2", "region2",
      "region1", "region3", "region3"
    ),
    age = c(1, 2, 3, 4, 5, 6, 7)
  )
  test_combined_gt <- cbind(test_gt, new_df)
  expect_true(inherits(test_combined_gt, "gen_tbl"))
  # check that the loci info is there
  expect_equal(nrow(show_loci(test_combined_gt)), nrow(test_loci))
})
