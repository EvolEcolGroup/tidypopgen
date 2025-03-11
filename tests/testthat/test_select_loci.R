test_that("select_loci subsets correctly", {
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c"),
    population = c("pop1", "pop1", "pop2")
  )
  test_genotypes <- rbind(
    c(1, 1, 0, 1, 1, 2),
    c(2, 1, 0, NA, 0, NA),
    c(2, 2, 0, 0, 1, NA)
  )
  test_loci <- data.frame(
    name = c(paste0("rs", 1:4), paste0("x", 1:2)),
    chromosome = c(1, 1, 1, 1, 2, 2),
    position = c(3, 5, 65, 343, 23, 456),
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

  # select snps with an rs
  test_gt_sub <- test_gt %>% select_loci(starts_with("rs"))
  expect_true(!any(c("x1", "x2") %in% loci_names(test_gt_sub)))
  # subsetting by id with reordering
  test_gt_sub <- test_gt %>% select_loci(c(3, 1, 5))
  expect_identical(c("rs3", "rs1", "x1"), loci_names(test_gt_sub))
  # get everything
  test_gt_sub <- test_gt %>% select_loci(everything())
  expect_identical(test_gt, test_gt_sub)
  # use 2:4 range expressions
  test_gt_sub <- test_gt %>% select_loci(2:4)
  expect_identical(test_loci$name[2:4], loci_names(test_gt_sub))
})
