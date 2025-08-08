skip_if_not_installed("hierfstat")

test_that("loci_pi computes correctly", {
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
  pi <- function(x) {
    n <- colSums(!is.na(x)) * 2
    c_0 <- colSums(x, na.rm = TRUE)
    c_1 <- n - c_0
    pi_est <- c_0 * c_1 / (n * (n - 1) / 2)
    return(pi_est)
  }
  pi_by_hand <- pi(test_genotypes)
  pi_gt <- loci_pi(test_gt$genotypes)
  expect_true(all(pi_by_hand == pi_gt))

  # compare to hierfstat (we don't give L, so it is just the sum of loci pi)
  pi_hierfstat <- hierfstat::pi.dosage(test_genotypes)
  expect_equal(sum(pi_gt), pi_hierfstat)


  # repeat the tests for a subset of the data
  # remove the 2nd individual and the 3rd and 5th snp
  test_genotypes_subset1 <- test_genotypes[-2, c(-3, -5)]
  test_gt_subset1 <- test_gt %>%
    filter(id != "b") %>%
    select_loci(c(-3, -5))
  pi_by_hand <- pi(test_genotypes_subset1)
  pi_gt <- loci_pi(test_gt_subset1$genotypes)
  expect_true(all(pi_by_hand == pi_gt))

  # now repeat with multiple blocks of snps
  pi_chunked <- loci_pi(test_gt_subset1, block_size = 2)
  expect_true(all(pi_gt == pi_chunked))

  # return NA if we have a single individual
  test_gt_single <- test_gt %>% filter(id == "a")
  pi_single <- test_gt_single %>% loci_pi()
  expect_true(all(is.na(pi_single)))
  # length should be the same as the number of loci
  expect_equal(length(pi_single), nrow(test_loci))
})

test_that("loci_pi on grouped tibbles", {
  test_genotypes <- rbind(
    c(1, 1, 0, 1, 1, 0),
    c(2, 1, 0, NA, 0, 0),
    c(2, NA, 0, 0, 1, 1),
    c(1, 0, 0, 1, 0, 0),
    c(1, 2, 0, 1, 2, 1),
    c(0, 0, 0, 0, NA, 1),
    c(0, 1, 1, 0, 1, NA)
  )
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c", "d", "e", "f", "g"),
    population = c("pop1", "pop1", "pop2", "pop2", "pop1", "pop3", "pop3")
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
  test_gt <- test_gt %>% group_by(population)

  # compute using .grouped_df method
  list <- loci_pi(test_gt, type = "list")
  matrix <- loci_pi(test_gt, type = "matrix")
  tidy <- loci_pi(test_gt, type = "tidy")
  expect_equal(list[1][[1]], as.vector(matrix[, "pop1"]))
  expect_equal(rownames(matrix), show_loci(test_gt)$name)
  expect_equal(colnames(matrix), group_keys(test_gt)$population)
  tidy_pop1 <- tidy %>%
    filter(group == "pop1") %>%
    select(value)
  expect_equal(list[1][[1]], tidy_pop1$value)

  # subset
  test_gt_subset <- test_gt %>% select_loci(c(1, 2, 3, 4))
  list <- loci_pi(test_gt_subset, type = "list")
  matrix <- loci_pi(test_gt_subset, type = "matrix")
  tidy <- loci_pi(test_gt_subset, type = "tidy")
  expect_equal(list[1][[1]], as.vector(matrix[, "pop1"]))
  expect_equal(rownames(matrix), show_loci(test_gt_subset)$name)
  expect_equal(colnames(matrix), group_keys(test_gt_subset)$population)
  tidy_pop1 <- tidy %>%
    filter(group == "pop1") %>%
    select(value)
  expect_equal(list[1][[1]], tidy_pop1$value)

  # compute by using group map
  pi_map <- test_gt %>% group_map(.f = ~ loci_pi(.x))
  # use fast cpp code (limit cores to 2)
  pi_grp <- test_gt %>% loci_pi(n_cores = 2)
  all.equal(pi_map, pi_grp)

  # and now with reframe
  loci_pi_reframe <-
    test_gt %>% reframe(missing = loci_pi(genotypes))
  loci_pi_direct <- test_gt %>%
    loci_pi(n_cores = 2) %>%
    arrange(group)
  expect_equal(loci_pi_reframe$missing, loci_pi_direct$value)
  # check that the direct method can take a column genotypes
  loci_pi_direct2 <- test_gt %>%
    loci_pi(genotypes) %>%
    arrange(group)
  expect_equal(loci_pi_reframe$missing, loci_pi_direct2$value)

  # now repeat with multiple blocks of snps
  pi_grp_chunked <- test_gt %>%
    loci_pi(n_cores = 2, block_size = 2)
  expect_true(all.equal(pi_grp, pi_grp_chunked))

  # return NA if we have a single individual
  test_gt_single <- test_gt %>% filter(id == "a")
  pi_single <- test_gt_single %>% loci_pi()
  expect_true(all(is.na(pi_single)))
  # length should be the same as the number of loci
  expect_equal(length(pi_single), nrow(test_loci))

  # test a second grouping variable
  test_gt$region <- c("a", "a", "b", "b", "a", "b", "b")
  test_gt <- test_gt %>% group_by(population, region)
  expect_error(
    test_gt %>% loci_pi(),
    "only works with one grouping variable"
  )
})
