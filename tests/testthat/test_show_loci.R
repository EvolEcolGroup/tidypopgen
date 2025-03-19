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

  # chromosome will be converted to character by gen_tibble
  test_loci$chromosome <- as.character(test_loci$chromosome)
  # check that we retrieve the info we put in (as a tibble)
  expect_equal(
    show_loci(test_gt$genotypes) %>% select(c(-big_index, -chr_int)),
    as_tibble(test_loci)
  )
  # now change it directly on the genotype column
  test_loci2 <- test_loci %>% dplyr::mutate(chromosome = "new")
  show_loci(test_gt$genotypes) <- test_loci2
  expect_identical(show_loci(test_gt), as_tibble(test_loci2))
  test_loci3 <- test_loci %>% dplyr::mutate(chromosome = "newer")
  show_loci(test_gt) <- test_loci3
  expect_identical(show_loci(test_gt), as_tibble(test_loci3))
  # with some proper dplyr
  show_loci(test_gt) <- show_loci(test_gt) %>% mutate(chromosome = "old")
  expect_true(all(show_loci(test_gt)$chromosome == "old"))
  test_loci3 <- test_loci3[-1, ]
  expect_error(
    show_loci(test_gt) <- test_loci3,
    "the replacement loci tibble does"
  )

  # try changing show_loci()$chromosome to integer
  test_loci4 <- test_loci %>%
    dplyr::mutate(chromosome = as.integer(c(1, 1, 2, 2, 3, 3)))
  show_loci(test_gt$genotypes) <- test_loci4
  # check method corrects it to a character
  expect_true(is.character(show_loci(test_gt)$chromosome))
  test_loci5 <- test_loci %>%
    dplyr::mutate(chromosome = as.integer(c(1, 2, 3, 4, 5, 6)))
  show_loci(test_gt) <- test_loci5
  # check method corrects it to a character
  expect_true(is.character(show_loci(test_gt)$chromosome))
})
