test_that("show_loci gets and sets information", {
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

  # change the chromosome names to a character with "chr" prefix
  show_loci(test_gt)$chromosome <- paste0("chr", rep(1:2, each = 3))
  # `show_loci()<-` will reset chromosome to factor; prefix "chr" is preserved
  expect_true(is.factor(show_loci(test_gt)$chromosome))
  expect_equal(levels(show_loci(test_gt)$chromosome), c("chr1", "chr2"))
  expect_equal(
    show_loci(test_gt)$chromosome,
    as.factor(c("chr1", "chr1", "chr1", "chr2", "chr2", "chr2"))
  )

  # `show_loci()<-` fails if replacement tibble has too few rows
  test_loci2 <- test_loci[-1, ]
  expect_error(
    show_loci(test_gt) <- test_loci2,
    "the replacement loci tibble does"
  )

  # `show_loci()<-` fails if replacement tibble has incomplete columns
  test_loci3 <- test_loci[, -1]
  expect_error(
    show_loci(test_gt) <- test_loci3,
    "loci must have the following columns"
  )

  # try changing show_loci()$chromosome to an integer
  show_loci(test_gt)$chromosome <- as.integer(c(1, 1, 2, 2, 3, 3))
  # check method corrects it to a factor
  expect_true(is.factor(show_loci(test_gt)$chromosome))
  expect_equal(levels(show_loci(test_gt)$chromosome), c("1", "2", "3"))
  expect_equal(
    show_loci(test_gt)$chromosome,
    as.factor(c("1", "1", "2", "2", "3", "3"))
  )
})
