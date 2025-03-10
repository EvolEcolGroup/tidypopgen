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


test_that("indiv_missingness computes correctly", {
  sum_na <- function(x) {
    sum(is.na(x))
  }
  # feeding the genotypes directly
  expect_true(all(
    indiv_missingness(test_gt$genotypes) ==
      apply(test_genotypes, 1, sum_na) / ncol(test_genotypes)
  ))
  # passing tibble
  expect_true(all(
    indiv_missingness(test_gt) ==
      apply(test_genotypes, 1, sum_na) / ncol(test_genotypes)
  ))
  # now using block_size to chunk the operation
  expect_true(all(
    indiv_missingness(test_gt, block_size = 2) ==
      apply(test_genotypes, 1, sum_na) / ncol(test_genotypes)
  ))
})


test_that("as_counts switch", {
  qc_rates <- indiv_missingness(test_gt, as_counts = FALSE)
  qc_counts <- indiv_missingness(test_gt, as_counts = TRUE)

  # All rates should be less than 1 - a proportion
  expect_true(all(qc_rates < 1))

  # Counts should be 0,2,1
  expect_equal(qc_counts, c(0, 2, 1))

  # They should not be identical outputs
  expect_false(identical(qc_rates, qc_counts))
})

test_genotypes_more_na <- rbind(
  c(1, 1, 0, 1, 1, 2),
  c(2, NA, NA, NA, NA, NA),
  c(2, 2, 0, 0, 1, NA)
)

test_gt2 <- gen_tibble(
  x = test_genotypes_more_na,
  loci = test_loci,
  indiv_meta = test_indiv_meta,
  quiet = TRUE
)


test_that("indiv_missingness on a subset", {
  # create a subset gt
  test_gt2_subset <- test_gt2 %>% select_loci(c(2, 3, 4, 5, 6))
  test_genotypes_more_na_subset <- test_genotypes_more_na[, c(2:6)]

  sum_na <- function(x) {
    sum(is.na(x))
  }
  # feeding the genotypes directly
  expect_true(all(
    indiv_missingness(test_gt2_subset$genotypes) ==
      round(
        apply(test_genotypes_more_na_subset, MARGIN = 1, sum_na) /
          ncol(test_genotypes_more_na_subset),
        digits = 1
      )
  ))
})
