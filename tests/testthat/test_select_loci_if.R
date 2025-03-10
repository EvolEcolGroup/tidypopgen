test_that("select_loci_if subsets correctly", {
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

  # raw frequencies
  freq <- test_gt %>% loci_maf()
  # remove monomorphic
  test_gen_sub <- test_gt %>% select_loci_if(loci_maf(genotypes) != 0)
  expect_true(!any(c("rs3", "x2") %in% show_loci(test_gen_sub)$name))
  # same subsetting with .data
  expect_identical(
    test_gen_sub,
    test_gt %>% select_loci_if(loci_maf(.data$genotypes) != 0)
  )
  # same manually
  criterion <- c(TRUE, TRUE, FALSE, TRUE, TRUE, FALSE)
  expect_identical(
    test_gen_sub,
    test_gt %>% select_loci_if(criterion)
  )
  # treat NA as if it was FALSE
  criterion_na <- c(TRUE, TRUE, NA, TRUE, TRUE, FALSE)
  expect_identical(
    test_gen_sub,
    test_gt %>% select_loci_if(criterion_na)
  )

  # Test error warnings: boolean criteria too short
  criterion_na <- c(TRUE, TRUE, NA, TRUE, TRUE)
  expect_error(
    test_gt %>% select_loci_if(criterion_na),
    ".sel_logical should be the same length "
  )

  # Test error warnings: selection criteria not boolean
  rsID <- c("rs1", "rs2", "rs3", "rs4", "x1", "x2") # nolint
  expect_error(
    test_gt %>% select_loci_if(rsID),
    ".sel_logical should be a logical boolean vector"
  )
})
