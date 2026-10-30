skip_if_not_installed("hierfstat")

# reference implementation of Tajima's D with missing data, computed directly
# from the dosage matrix. Following Bailey et al. (2025) and pixy, each locus
# uses its own sample size (the number of valid allele copies) to estimate pi
# and Watterson's theta; the standard deviation of Tajima (1989) is computed
# for the (rounded) mean sample size across loci. hierfstat::TajimaD.dosage()
# uses the same sample size for all loci, so it can only be used to validate
# the case with no missing data. See also test_pixy_tajimas_d.R for a direct
# comparison with pixy.
ref_tajimas_d <- function(genotypes, ploidy) {
  # sum_{i=1}^{n-1} 1/i^power
  a_n <- function(n, power = 1) sum(1 / (seq_len(n - 1))^power)
  # number of valid allele copies, and counts of each allele, at each locus
  n_x <- apply(genotypes, 2, function(col) sum(ploidy[!is.na(col)]))
  c_1 <- colSums(genotypes, na.rm = TRUE)
  c_0 <- n_x - c_1
  # only use loci with at least one valid allele copy
  observed <- n_x > 0
  n_x <- n_x[observed]
  c_0 <- c_0[observed]
  c_1 <- c_1[observed]
  # pi at each locus (0 if there are fewer than 2 allele copies)
  pi_x <- ifelse(n_x > 1, c_0 * c_1 / (n_x * (n_x - 1) / 2), 0)
  # segregating sites, and watterson's theta weighting each segregating site
  # by its own sample size
  is_seg <- c_0 > 0 & c_1 > 0
  seg <- sum(is_seg)
  theta_w <- sum(1 / sapply(n_x[is_seg], a_n))
  # standard deviation for the mean sample size
  n <- round(mean(n_x))
  a1 <- a_n(n)
  a2 <- a_n(n, power = 2)
  e1 <- ((n + 1) / (3 * (n - 1)) - 1 / a1) / a1
  e2_num <- 2 * (n^2 + n + 3) / (9 * n * (n - 1)) - (n + 2) / (n * a1) +
    a2 / (a1^2)
  e2 <- e2_num / (a1^2 + a2)
  (sum(pi_x) - theta_w) / sqrt(e1 * seg + e2 * seg * (seg - 1))
}

test_that("pop_tajimas_d computes correctly", {
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

  taj_gt <- pop_tajimas_d(test_gt$genotypes)

  # with missing data, each locus uses its own sample size
  expect_equal(taj_gt, ref_tajimas_d(test_genotypes, rep(2, 3)))
  # compare to hierfstat on the loci with no missing data
  # (we don't give L, so it is just the sum of loci pi)
  taj_gt_complete <- pop_tajimas_d(select_loci(test_gt, c(-4, -6))$genotypes)
  taj_hierfstat <- hierfstat::TajimaD.dosage(test_genotypes[, c(-4, -6)])
  expect_equal(taj_gt_complete, taj_hierfstat)


  # repeat the tests for a subset of the data
  # remove the 2nd individual and the 3rd and 5th snp
  test_genotypes_subset1 <- test_genotypes[-2, c(-3, -5)]
  test_gt_subset1 <- test_gt %>%
    filter(id != "b") %>%
    select_loci(c(-3, -5))
  taj_gt <- pop_tajimas_d(test_gt_subset1$genotypes)
  expect_equal(taj_gt, ref_tajimas_d(test_genotypes_subset1, rep(2, 2)))
  # the last locus has missing data, so we remove it to compare to hierfstat
  taj_hierfstat <- hierfstat::TajimaD.dosage(test_genotypes_subset1[, -4])
  expect_equal(
    pop_tajimas_d(select_loci(test_gt_subset1, -4)$genotypes),
    taj_hierfstat
  )
  # now repeat with multiple blocks of snps
  taj_chunked <- pop_tajimas_d(test_gt_subset1, block_size = 2)
  expect_true(taj_gt == taj_chunked)
})

test_that("pop_tajimas_d on grouped tibbles", {
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
  # compute by using group map
  taj_map <- test_gt %>% group_map(.f = ~ pop_tajimas_d(.x))
  # use fast cpp code (limit cores to 2)
  taj_grp <- test_gt %>% pop_tajimas_d(n_cores = 2)
  expect_true(all.equal(taj_map, taj_grp))
  # now repeat with multiple blocks of snps
  taj_grp_chunked <- test_gt %>%
    pop_tajimas_d(n_cores = 2, block_size = 2)
  expect_true(all.equal(taj_grp, taj_grp_chunked))
})

test_that("pop_tajimas_d ignores loci with missing pi", {
  test_genotypes <- rbind(
    c(1, 1, 0, 1, NA, 0),
    c(2, 1, 0, NA, NA, 0),
    c(2, NA, 0, 0, 1, 1),
    c(1, 0, 0, 1, 0, 0),
    c(1, 2, 0, 1, NA, 1),
    c(0, 0, 0, 0, NA, 1),
    c(0, 1, 1, 0, NA, NA)
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
  # the 5th locus is missing in pop1 and pop3, so pi is NA for those groups
  pi_by_grp <- loci_pi(test_gt, type = "list")
  expect_true(is.na(pi_by_grp[[1]][5]))
  expect_true(is.na(pi_by_grp[[3]][5]))
  # tajima's d is computed ignoring that locus
  taj_grp <- test_gt %>% pop_tajimas_d(n_cores = 2)
  expect_false(any(is.na(unlist(taj_grp))))
  taj_grp_subset <- test_gt %>%
    select_loci(-5) %>%
    pop_tajimas_d(n_cores = 2)
  expect_equal(taj_grp[c(1, 3)], taj_grp_subset[c(1, 3)])
  # pop2 has data for the 5th locus, so it is used
  expect_false(taj_grp[[2]] == taj_grp_subset[[2]])
  # same for ungrouped data
  test_gt_pop1 <- test_gt %>%
    filter(population == "pop1") %>%
    ungroup()
  expect_equal(
    pop_tajimas_d(test_gt_pop1),
    pop_tajimas_d(test_gt_pop1 %>% select_loci(-5))
  )
})

test_that("pop_tajimas_d counts pseudohaploid allele copies", {
  pseudohaploid_genotypes <- rbind(
    c(0, 2, NA),
    c(1, 2, 0),
    c(0, NA, 2),
    c(1, 0, 0)
  )
  pseudohaploid_meta <- data.frame(
    id = letters[1:4],
    population = c("popA", "popA", "popB", "popB")
  )
  pseudohaploid_loci <- data.frame(
    name = paste0("rs", 1:3),
    chromosome = 1,
    position = 1:3,
    genetic_dist = 0,
    allele_ref = "A",
    allele_alt = "T"
  )
  pseudohaploid_gt <- gen_tibble(
    x = pseudohaploid_genotypes,
    loci = pseudohaploid_loci,
    indiv_meta = pseudohaploid_meta,
    quiet = TRUE
  ) %>%
    gt_pseudohaploid(test_n_loci = NULL) %>%
    dplyr::group_by(population)

  taj_grp <- pop_tajimas_d(pseudohaploid_gt)
  pi_by_group <- loci_pi(pseudohaploid_gt, type = "list")
  # number of valid allele copies at each locus for each group (pseudohaploids
  # contribute one copy, diploids two)
  ploidy <- indiv_ploidy(pseudohaploid_gt)
  n_by_group <- lapply(
    split(seq_along(ploidy), pseudohaploid_meta$population),
    function(rows) {
      valid <- !is.na(pseudohaploid_genotypes[rows, , drop = FALSE])
      colSums(valid * ploidy[rows])
    }
  )

  expected <- Map(tajimas_d_from_pi_vec, pi_by_group, n_by_group)
  expect_equal(unlist(taj_grp), unname(unlist(expected)))
  # same as computing on each group, and as a window covering all loci
  taj_map <- pseudohaploid_gt %>% group_map(.f = ~ pop_tajimas_d(.x))
  expect_equal(taj_map, taj_grp)
  window_taj <- windows_pop_tajimas_d(
    pseudohaploid_gt,
    type = "list",
    window_size = 3,
    step_size = 3
  )
  expect_equal(
    unname(vapply(window_taj, function(x) x$stat[1], numeric(1))),
    unlist(taj_grp)
  )
})

test_that("pop_tajimas_d works on polyploid data", {
  test_indiv_meta <- data.frame(
    id = letters[1:6],
    population = c("popA", "popA", "popA", "popB", "popB", "popB")
  )
  test_genotypes <- rbind(
    c(1, 4, 0, 2),
    c(4, 3, 0, 0),
    c(0, 2, 1, 4),
    c(2, 0, 4, NA),
    c(1, 1, 2, 2),
    c(3, 0, 4, 1)
  )
  test_ploidy <- rep(4, 6)
  test_loci <- data.frame(
    name = paste0("rs", 1:4),
    chromosome = c(1, 1, 2, 2),
    position = c(3, 5, 23, 456),
    genetic_dist = as.double(rep(0, 4)),
    allele_ref = c("A", "T", "C", "G"),
    allele_alt = c("T", "C", "G", "A")
  )
  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    ploidy = test_ploidy,
    quiet = TRUE
  )

  # ungrouped
  taj_gt <- pop_tajimas_d(test_gt)
  expect_equal(taj_gt, ref_tajimas_d(test_genotypes, test_ploidy))
  # grouped
  test_gt <- test_gt %>% group_by(population)
  taj_map <- test_gt %>% group_map(.f = ~ pop_tajimas_d(.x))
  taj_grp <- test_gt %>% pop_tajimas_d(n_cores = 2)
  expect_equal(taj_map, taj_grp)
  expect_equal(taj_grp[[1]], ref_tajimas_d(test_genotypes[1:3, ], rep(4, 3)))
  expect_equal(taj_grp[[2]], ref_tajimas_d(test_genotypes[4:6, ], rep(4, 3)))
})

test_that("pop_tajimas_d uses the sample size of each locus", {
  test_genotypes <- rbind(
    c(1, 1, NA, 1, 0),
    c(2, 1, 1, NA, 0),
    c(0, 2, NA, 0, 1),
    c(1, 0, NA, 1, NA)
  )
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c", "d"),
    population = c("pop1", "pop1", "pop1", "pop1")
  )
  test_loci <- data.frame(
    name = paste0("rs", 1:5),
    chromosome = c(1, 1, 1, 2, 2),
    position = c(3, 5, 65, 23, 456),
    genetic_dist = as.double(rep(0, 5)),
    allele_ref = c("A", "T", "C", "G", "C"),
    allele_alt = c("T", "C", "G", "C", "G")
  )
  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  )
  # compare to the reference implementation
  expect_equal(pop_tajimas_d(test_gt), ref_tajimas_d(test_genotypes, rep(2, 4)))
  # the 3rd locus is typed in a single heterozygous individual (pi = 1), and
  # it should be counted as segregating
  expect_equal(loci_pi(test_gt)[3], 1)
  expect_equal(unname(tajimas_d_locus_terms(1, 2)[, "seg"]), 1)
  # with no missing data, all loci have the same sample size, and we recover
  # the classic Tajima's D with a fixed sample size
  test_gt_complete <- test_gt %>% select_loci(c(1, 2))
  expect_equal(
    pop_tajimas_d(test_gt_complete),
    tajimas_d_from_pi_vec(loci_pi(test_gt_complete), 8)
  )
  expect_equal(
    pop_tajimas_d(test_gt_complete),
    hierfstat::TajimaD.dosage(test_genotypes[, c(1, 2)])
  )
})
