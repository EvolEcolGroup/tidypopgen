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

test_gt <- gen_tibble(x = test_genotypes, loci = test_loci, indiv_meta = test_indiv_meta, quiet = TRUE)


test_that("loci_missingness", {
  # na counts
  count_na <- function(x) {
    sum(is.na(x))
  }
  n_na <- apply(test_genotypes, 2, count_na)
  expect_true(all(loci_missingness(test_gt$genotypes, as_counts = TRUE) == n_na))
  # convert to frequencies
  expect_true(all(loci_missingness(test_gt$genotypes) == n_na / nrow(test_genotypes)))

  # create a subset gt
  test_gt_subset1 <- test_gt %>% filter(id != "a")
  # subset genotypes
  test_genotypes_subset1 <- test_genotypes[-1, ]

  n_na <- apply(test_genotypes_subset1, 2, count_na)
  expect_true(all(loci_missingness(test_gt_subset1$genotypes, as_counts = TRUE) == n_na))
  # convert to frequencies
  expect_true(all(loci_missingness(test_gt_subset1$genotypes) == n_na / nrow(test_genotypes_subset1)))
})

test_that("loci_missingness computes correctly", {
  # na counts
  count_na <- function(x) {
    sum(is.na(x))
  }
  n_na <- apply(test_genotypes, 2, count_na)
  expect_true(all(loci_missingness(test_gt$genotypes, as_counts = TRUE) == n_na))
  # convert to frequencies
  expect_true(all(loci_missingness(test_gt$genotypes) == n_na / nrow(test_genotypes)))

  # now using block_size to chunk the operation
  expect_true(all(loci_missingness(test_gt, block_size = 2) ==
    n_na / nrow(test_genotypes)))
})

test_that("loci_missingness on grouped tibble", {
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

  test_gt <- gen_tibble(x = test_genotypes, loci = test_loci, indiv_meta = test_indiv_meta, quiet = TRUE)
  test_gt <- test_gt %>% group_by(population)
  # compute by using group map
  loci_miss_map <- test_gt %>% group_map(.f = ~ loci_missingness(.x))
  # use fast cpp code (limit cores to 2)
  loci_miss_grp <- test_gt %>% loci_missingness(n_cores = 2)
  expect_true(all.equal(loci_miss_map, loci_miss_grp))
})
