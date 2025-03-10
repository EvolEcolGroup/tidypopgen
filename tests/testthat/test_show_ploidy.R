test_that("ploidy works correctly", {
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
    chromosome = as.integer(c(1, 1, 1, 1, 2, 2)),
    position = as.integer(c(3, 5, 65, 343, 23, 456)),
    genetic_dist = as.integer(rep(0, 6)),
    allele_ref = c("A", "T", "C", "G", "C", "T"),
    allele_alt = c("T", "C", NA, "C", "G", "A")
  )
  test_gt <- gen_tibble(
    test_genotypes,
    indiv_meta = test_indiv_meta,
    loci = test_loci,
    backingfile = tempfile(),
    quiet = TRUE
  )
  # we expect a diploid by default
  expect_true(show_ploidy(test_gt) == 2)
  expect_true(show_ploidy(test_gt$genotypes) == 2)
  # get ploidy for each individual
  expect_true(all(indiv_ploidy(test_gt) == c(2, 2, 2)))
  expect_true(all(indiv_ploidy(test_gt$genotypes) == c(2, 2, 2)))
  # now let's create a tetraploid for individual 2
  test_genotypes <- rbind(
    c(1, 1, 0, 1, 1, 2),
    c(4, 1, 3, NA, 0, NA),
    c(2, 2, 0, 0, 1, NA)
  )
  # throw an error if we request mixed ploidy but don't have the ploidy info
  expect_error(
    gen_tibble(
      test_genotypes,
      indiv_meta = test_indiv_meta,
      loci = test_loci,
      backingfile = tempfile(),
      ploidy = 0,
      quiet = TRUE
    ),
    "max ploidy "
  )
  # add ploidy info for each individual
  test_gt <- gen_tibble(
    test_genotypes,
    indiv_meta = test_indiv_meta,
    loci = test_loci,
    backingfile = tempfile(),
    ploidy = c(2, 4, 2),
    quiet = TRUE
  )
  # check that we have the correct genotypes
  expect_true(all.equal(
    test_genotypes,
    show_genotypes(test_gt)
  ))
  expect_true(show_ploidy(test_gt) == 0)
  # get ploidy for each individual
  expect_true(all(indiv_ploidy(test_gt) == c(2, 4, 2)))
  expect_true(all(indiv_ploidy(test_gt$genotypes) == c(2, 4, 2)))

  # expect warning when trying to do something with a mixed ploidy object
  expect_warning(test_gt %>% loci_maf(), "this function ")

  # missingness
  count_na <- function(x) {
    sum(is.na(x))
  }
  freq_na <- apply(test_genotypes, 2, count_na) / nrow(test_genotypes)
  expect_true(all(loci_missingness(test_gt$genotypes) == freq_na))
})
