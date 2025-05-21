test_genotypes <- rbind(
  c(1, 1, 0, 1, 1, 0),
  c(2, 1, 0, NA, 0, 0),
  c(2, NA, 0, 0, 1, 1),
  c(1, 0, 0, 1, 0, 0),
  c(1, 2, 0, 1, 2, 1),
  c(0, 0, 0, 0, NA, 1),
  c(0, 1, 1, 0, 1, NA),
  c(1, 1, 2, 1, 1, 1),
  c(1, 2, 1, 0, 1, 1),
  c(0, 1, 1, 1, NA, 2)
)
test_indiv_meta <- data.frame(
  id = c("a", "b", "c", "d", "e", "f", "g", "h", "i", "k"),
  population = c(
    "pop1", "pop1", "pop2", "pop2", "pop1", "pop3", "pop3",
    "pop4", "pop4", "pop4"
  )
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

test_that("requires a grouped gen_tibble", {
  expect_error(
    test_gt %>% nwise_pop_pbs(fst_method = "Hudson"),
    ".x should be a grouped gen_tibble"
  )
})
test_that("error with multiple grouping variables", {
  # test a second grouping variable
  test_gt$region <- c("a", "a", "b", "b", "a", "b", "b", "a", "b", "a")
  test_gt <- test_gt %>% group_by(population, region)
  expect_error(
    test_gt %>% nwise_pop_pbs(fst_method = "Hudson"),
    "only works with one grouping variable"
  )
})
test_that("nwise_pop_pbs works correctly", {
  test_pbs <- test_gt %>%
    group_by(population) %>%
    nwise_pop_pbs(fst_method = "Hudson", type = "matrix")
  # expect 24 columns
  expect_equal(ncol(test_pbs), 24)
  # @TODO we need some more meaningul tests here
  # maybe check sk.allele if it does pbs
})

test_that("type argument delivers correct objects", {
  test_pbs_matrix <- test_gt %>%
    group_by(population) %>%
    nwise_pop_pbs(fst_method = "Hudson", type = "matrix")
  expect_true(is.matrix(test_pbs_matrix))

  test_pbs_tidy <- test_gt %>%
    group_by(population) %>%
    nwise_pop_pbs(fst_method = "Hudson", type = "tidy")
  expect_true(is.data.frame(test_pbs_tidy))

  # Compare
  pop1_pop2_pop3_tidy <-
    subset(test_pbs_tidy, test_pbs_tidy$stat_name == "pbs_pop1.pop2.pop3")
  expect_equal(pop1_pop2_pop3_tidy$value, as.vector(test_pbs_matrix[, 1]))
})
