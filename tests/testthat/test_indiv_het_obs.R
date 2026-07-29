# limit number of threads for tests
data.table::setDTthreads(2)
if (rlang::is_installed("RhpcBLASctl")) {
  RhpcBLASctl::blas_set_num_threads(2)
  RhpcBLASctl::omp_set_num_threads(2)
}

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

test_that("indiv_het_obs computes correctly for tetraploid data", {
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c"),
    population = c("pop1", "pop1", "pop2")
  )
  # dosages of the alt allele, range 0-4 for a tetraploid
  test_genotypes_tetra <- rbind(
    c(1, 4, 0, 2, 3, 4),
    c(4, 4, 0, 0, 1, NA),
    c(0, 2, 0, 4, 4, NA)
  )
  test_loci <- data.frame(
    name = paste0("rs", 1:6),
    chromosome = c(1, 1, 1, 1, 2, 2),
    position = c(3, 5, 65, 343, 23, 456),
    genetic_dist = as.double(rep(0, 6)),
    allele_ref = c("A", "T", "C", "G", "C", "T"),
    allele_alt = c("T", "C", NA, "C", "G", "A")
  )

  test_gt_tetra <- gen_tibble(
    x = test_genotypes_tetra,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    ploidy = 4,
    quiet = TRUE
  )
  expect_true(all(indiv_ploidy(test_gt_tetra) == 4))

  # heterozygous whenever the dosage is neither 0 nor the ploidy (4)
  expected_het <- rowMeans(
    test_genotypes_tetra > 0 & test_genotypes_tetra < 4,
    na.rm = TRUE
  )
  expect_equal(indiv_het_obs(test_gt_tetra), expected_het)
  expect_equal(indiv_het_obs(test_gt_tetra$genotypes), expected_het)

  # counts
  counts_mat <- as.matrix(
    data.frame(
      het_n = apply(
        test_genotypes_tetra > 0 & test_genotypes_tetra < 4,
        1, sum,
        na.rm = TRUE
      ),
      na_n = apply(is.na(test_genotypes_tetra), 1, sum, na.rm = TRUE)
    )
  )
  expect_identical(
    indiv_het_obs(test_gt_tetra$genotypes, as_counts = TRUE),
    counts_mat
  )
})

test_that("indiv_het_obs computes correctly for mixed-ploidy data", {
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c", "d"),
    population = c("pop1", "pop1", "pop2", "pop2")
  )
  # individuals a,b are diploid; c,d are tetraploid
  test_genotypes_mixed <- rbind(
    c(1, 2, 0, 1, NA, 2),
    c(0, 2, 1, 0, 2, 1),
    c(4, 4, 0, 2, 3, NA),
    c(0, 1, 4, 4, 4, 2)
  )
  test_ploidy <- c(2, 2, 4, 4)
  test_loci <- data.frame(
    name = paste0("rs", 1:6),
    chromosome = c(1, 1, 1, 1, 2, 2),
    position = c(3, 5, 65, 343, 23, 456),
    genetic_dist = as.double(rep(0, 6)),
    allele_ref = c("A", "T", "C", "G", "C", "T"),
    allele_alt = c("T", "C", NA, "C", "G", "A")
  )

  test_gt_mixed <- gen_tibble(
    x = test_genotypes_mixed,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    ploidy = test_ploidy,
    quiet = TRUE
  )
  expect_identical(indiv_ploidy(test_gt_mixed), test_ploidy)
  # mixed ploidy at the gen_tibble level
  expect_true(show_ploidy(test_gt_mixed) == 0)

  # each individual is evaluated against its own ploidy
  expected_het <- vapply(
    seq_len(nrow(test_genotypes_mixed)),
    function(i) {
      mean(
        test_genotypes_mixed[i, ] > 0 &
          test_genotypes_mixed[i, ] < test_ploidy[i],
        na.rm = TRUE
      )
    },
    numeric(1)
  )
  expect_equal(indiv_het_obs(test_gt_mixed), expected_het)
})
