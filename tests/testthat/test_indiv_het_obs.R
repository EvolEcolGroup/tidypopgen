test_that("indiv_het_obs computes correctly", {
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
    name = paste0("rs", 1:6),
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

  # feeding the list of SNPbin directly
  expect_true(all(
    indiv_het_obs(test_gt$genotypes) ==
      rowMeans(test_genotypes == 1, na.rm = TRUE)
  ))
  # passing tibble
  expect_true(all(
    indiv_het_obs(test_gt) == rowMeans(test_genotypes == 1, na.rm = TRUE)
  ))

  # get counts
  counts_mat <- as.matrix(
    data.frame(
      het_n = apply(test_genotypes == 1, 1, sum, na.rm = TRUE),
      na_n = apply(is.na(test_genotypes), 1, sum, na.rm = TRUE)
    )
  )
  expect_identical(
    indiv_het_obs(test_gt$genotypes, as_counts = TRUE),
    counts_mat
  )
})

test_that("indiv_het_obs returns 0's when all genotypes are homozygous", {
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c"),
    population = c("pop1", "pop1", "pop2")
  )
  test_genotypes_homozyg <- rbind(
    c(2, 2, 0, 0, 2, 0),
    c(2, 2, 0, 0, 2, 0),
    c(2, 2, 0, 0, 2, 0)
  )
  test_loci <- data.frame(
    name = paste0("rs", 1:6),
    chromosome = c(1, 1, 1, 1, 2, 2),
    position = c(3, 5, 65, 343, 23, 456),
    genetic_dist = as.double(rep(0, 6)),
    allele_ref = c("A", "T", "C", "G", "C", "T"),
    allele_alt = c("T", "C", NA, "C", "G", "A")
  )

  test_gt_homozyg <- gen_tibble(
    x = test_genotypes_homozyg,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  )

  # feeding the list of SNPbin directly
  expect_true(all(
    indiv_het_obs(test_gt_homozyg$genotypes) ==
      rowMeans(test_genotypes_homozyg == 1, na.rm = TRUE)
  ))

  # passing tibble
  expect_true(all(
    indiv_het_obs(test_gt_homozyg) ==
      rowMeans(test_genotypes_homozyg == 1, na.rm = TRUE)
  ))
})
