# tests for pop_het_obs() on polyploid and mixed-ploidy gen_tibbles.
# hierfstat (used to validate the diploid case in test_pop_basic_stats.R)
# is diploid-only, so here we validate against an independent reference
# implementation computed directly from the raw dosage matrix.

# reference implementation: for each locus/population, the proportion of
# non-missing individuals whose dosage is neither 0 nor their own ploidy.
ref_pop_het_obs <- function(genotypes, ploidy, population) {
  pop_levels <- unique(population)
  het <- genotypes > 0 & genotypes != matrix(
    ploidy,
    nrow = nrow(genotypes), ncol = ncol(genotypes)
  )
  out <- sapply(pop_levels, function(pop) { # nolint
    rows <- population == pop
    colMeans(het[rows, , drop = FALSE], na.rm = TRUE)
  })
  colnames(out) <- pop_levels
  out
}

test_that("pop_het_obs computes correctly for a single tetraploid population", {
  test_indiv_meta <- data.frame(
    id = letters[1:4],
    population = rep("popA", 4)
  )
  test_genotypes <- rbind(
    c(1, 4, 0, 2),
    c(4, 4, 0, 0),
    c(0, 2, 0, 4),
    c(2, 0, 4, NA)
  )
  test_ploidy <- rep(4, 4)
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
  ) %>% dplyr::group_by(population)

  expected <- ref_pop_het_obs(
    test_genotypes,
    test_ploidy,
    test_indiv_meta$population
  )

  test_het_obs <- test_gt %>% pop_het_obs(by_locus = TRUE)
  expect_equal(unname(test_het_obs), unname(expected))

  test_het_obs_mean <- test_gt %>% pop_het_obs(by_locus = FALSE)
  expect_equal(test_het_obs_mean, colMeans(expected, na.rm = TRUE))
})

test_that("pop_het_obs computes correctly for two tetraploid populations", {
  test_indiv_meta <- data.frame(
    id = letters[1:6],
    population = c("popA", "popA", "popA", "popB", "popB", "popB")
  )
  test_genotypes <- rbind(
    c(1, 4, 0, 2),
    c(4, 4, 0, 0),
    c(0, 2, 0, 4),
    c(2, 0, 4, 1),
    c(0, 0, 4, 3),
    c(4, 2, 2, NA)
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
  ) %>% dplyr::group_by(population)

  expected <- ref_pop_het_obs(
    test_genotypes,
    test_ploidy,
    test_indiv_meta$population
  )

  test_het_obs <- test_gt %>% pop_het_obs(by_locus = TRUE)
  expect_equal(unname(test_het_obs), unname(expected))
  expect_identical(colnames(test_het_obs), c("popA", "popB"))

  # genome-wide global estimate for polyploid data is the unweighted mean
  # of the population estimates (no dependency on pop_global_stats())
  test_het_obs_global <- test_gt %>%
    pop_het_obs(by_locus = TRUE, include_global = TRUE)
  expect_equal(
    unname(test_het_obs_global[, "global"]),
    rowMeans(expected, na.rm = TRUE)
  )
})

test_that("pop_het_obs handles a population mixing diploid and tetraploid individuals", { # nolint
  test_indiv_meta <- data.frame(
    id = letters[1:4],
    population = rep("popA", 4)
  )
  # a,b diploid; c,d tetraploid, all in the same population
  test_genotypes <- rbind(
    c(1, 2, 0, 1),
    c(0, 2, 1, 0),
    c(4, 4, 0, 2),
    c(0, 1, 4, 4)
  )
  test_ploidy <- c(2, 2, 4, 4)
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
  ) %>% dplyr::group_by(population)

  # each individual's own ploidy is used, so Ho is simply the proportion of
  # heterozygous individuals, unaffected by their ploidy differing
  expected <- ref_pop_het_obs(
    test_genotypes,
    test_ploidy,
    test_indiv_meta$population
  )
  test_het_obs <- test_gt %>% pop_het_obs(by_locus = TRUE)
  expect_equal(unname(test_het_obs), unname(expected))
})

test_that("pop_het_obs on polyploid data matches diploid results for ploidy = 2", { # nolint
  # sanity check: the new ploidy-agnostic kernel used for pop_het_obs()
  # must exactly reproduce the diploid results already validated against
  # hierfstat in test_pop_basic_stats.R
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c", "d", "e", "f", "g"),
    population = c("pop1", "pop1", "pop2", "pop2", "pop1", "pop3", "pop3")
  )
  test_genotypes <- rbind(
    c(1, 1, 0, 1, 1, 0),
    c(2, 1, 0, NA, 0, 0),
    c(2, NA, 0, 0, 1, 1),
    c(1, 0, 0, 1, 0, 0),
    c(1, 2, 0, 1, 2, 1),
    c(0, 0, 0, 0, NA, 1),
    c(0, 1, 1, 0, 1, NA)
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
  ) %>% dplyr::group_by(population)

  expected <- ref_pop_het_obs(
    test_genotypes,
    rep(2, nrow(test_genotypes)),
    test_indiv_meta$population
  )
  test_het_obs <- test_gt %>% pop_het_obs(by_locus = TRUE)
  expect_equal(unname(test_het_obs), unname(expected))
})
