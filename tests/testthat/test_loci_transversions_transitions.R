test_that("find transitions and transversions", {
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
  transv_bool <- c(TRUE, FALSE, NA, TRUE, TRUE, TRUE)
  expect_true(all.equal(loci_transversions(test_gt), transv_bool))
  expect_true(all.equal(loci_transitions(test_gt), !transv_bool))
})

test_that("check warning message for different alleles", {
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
    allele_ref = c("a", "t", "c", "g", "c", "t"),
    allele_alt = c("t", "c", NA, "c", "g", "a")
  )

  # Create a gen_tibble with alleles that are not "A" "T" "C" "G"
  # using valid_alleles
  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE,
    valid_alleles = c("a", "t", "c", "g")
  )

  testthat::expect_error(
    loci_transversions(test_gt),
    "valid alleles are A T C G 0 . but "
  )
  testthat::expect_error(
    loci_transitions(test_gt),
    "valid alleles are A T C G 0 . but "
  )
})
