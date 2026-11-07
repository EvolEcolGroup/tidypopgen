# Check that the jackknife of each group (a list as returned by
# jackknife_pop_tajimas_d(type = "list")) matches the jackknife computed by
# recomputing tajima's d from the genotypes of the loci left when leaving out
# each block, with the variance computed for the sample size of all loci
expect_jack_matches_subsets <- function(test_gt, block_id, jack_list) {
  # tajima's d of a gen_tibble (as a single group), with sample size n_var
  # in the variance (by default, the rounded mean sample size of its loci)
  tajimas_d_gt <- function(.x, n_var = NULL) {
    term_sums <- .tajimas_d_group_sums(.tajimas_d_sums_by_group(.x), 1)
    tajimas_d_from_term_sums(term_sums, n_var = n_var)
  }
  for (i_grp in seq_along(jack_list)) {
    # the individuals of the group (the groups are in the order of jack_list)
    pop_gt <- test_gt %>%
      dplyr::filter(dplyr::cur_group_id() == i_grp) %>%
      dplyr::ungroup()
    testthat::expect_equal(
      as.character(unique(pop_gt$population)),
      names(jack_list)[i_grp]
    )
    taj_full <- pop_tajimas_d(pop_gt)
    testthat::expect_equal(tajimas_d_gt(pop_gt), taj_full)
    term_full <- .tajimas_d_group_sums(.tajimas_d_sums_by_group(pop_gt), 1)
    n_var <- tajimas_d_n_var(term_full[, "n"], term_full[, "n_loci"])
    taj_loo <- vapply(sort(unique(block_id)), function(i_block) {
      pop_gt %>%
        select_loci(which(block_id != i_block)) %>%
        tajimas_d_gt(n_var = n_var)
    }, numeric(1))
    m <- as.vector(table(block_id))
    expected <- .jackknife_weighted(taj_full, matrix(taj_loo), matrix(m))
    testthat::expect_equal(jack_list[[i_grp]], expected)
  }
}

test_that("jackknife_pop_tajimas_d leave-one-out matches subsets", {
  test_gt <- load_example_gt("gen_tbl_gen_dist")
  block_id <- .jackknife_block_ids(show_loci(test_gt),
    block_size = 1, size_unit = "cM"
  )
  jack_list <- jackknife_pop_tajimas_d(test_gt,
    type = "list", block_size = 1, size_unit = "cM"
  )
  # estimate is the genome-wide tajima's d
  expect_equal(
    unname(vapply(jack_list, function(x) x$estimate, numeric(1))),
    unlist(pop_tajimas_d(test_gt))
  )
  expect_jack_matches_subsets(test_gt, block_id, jack_list)
})

test_that("jackknife_pop_tajimas_d uses the sample size of all loci", {
  # popA has missing genotypes, so that the sample sizes (in allele copies)
  # of its 4 blocks of 3 loci are 20,20,20 | 20,20,20 | 16,16,16 | 18,16,16:
  # the mean sample size over all loci (18.17) rounds to 18, but leaving out
  # the third or fourth block rounds it to 19
  set.seed(42)
  genotypes <- matrix(sample(0:2, 14 * 12, replace = TRUE), nrow = 14)
  genotypes[1:2, 7:9] <- NA
  genotypes[1, 10] <- NA
  genotypes[1:2, 11:12] <- NA
  test_gt <- gen_tibble(
    x = genotypes,
    loci = data.frame(
      name = paste0("rs", 1:12),
      chromosome = 1,
      position = 1:12,
      genetic_dist = 0,
      allele_ref = "A",
      allele_alt = "T"
    ),
    indiv_meta = data.frame(
      id = paste0("ind", 1:14),
      population = rep(c("popA", "popB"), c(10, 4))
    ),
    quiet = TRUE
  ) %>%
    dplyr::group_by(population)
  block_id <- .jackknife_block_ids(show_loci(test_gt),
    block_size = 3, size_unit = "snp"
  )
  # check that leaving out some blocks changes the rounding in popA
  pop_a <- test_gt %>%
    dplyr::ungroup() %>%
    dplyr::filter(population == "popA")
  block_sums <- .tajimas_d_group_sums(
    .tajimas_d_sums_by_group(pop_a,
      locus_first = block_id - 1L,
      locus_last = block_id - 1L,
      n_windows = max(block_id)
    ),
    1
  )
  expect_equal(block_sums[, "n"], c(60, 60, 48, 50))
  loo_n <- (sum(block_sums[, "n"]) - block_sums[, "n"]) /
    (sum(block_sums[, "n_loci"]) - block_sums[, "n_loci"])
  expect_equal(round(loo_n), c(18, 18, 19, 19))

  jack_list <- jackknife_pop_tajimas_d(test_gt,
    type = "list", block_size = 3, size_unit = "snp"
  )
  expect_equal(
    unname(vapply(jack_list, function(x) x$estimate, numeric(1))),
    unlist(pop_tajimas_d(test_gt))
  )
  expect_jack_matches_subsets(test_gt, block_id, jack_list)
})

test_that("the jackknife of tajimas d depends on the total sample size", {
  # sums for 10 blocks of 100 loci, with mean sample sizes of 19 or 18 allele
  # copies: the mean sample size over all loci (18.501) rounds to 19, but
  # leaving out a block with a mean of 19 rounds it to 18
  n_loci <- rep(100, 10)
  n_sum <- c(1901, rep(1900, 4), rep(1800, 5))
  seg <- c(40, 42, 38, 45, 41, 39, 43, 37, 44, 40)
  term_blocks <- cbind(
    k_hat = seg * c(0.21, 0.25, 0.23, 0.2, 0.22, 0.24, 0.26, 0.21, 0.23, 0.22),
    seg = seg,
    theta_w = seg / harmonic_sum(19),
    n = n_sum,
    n_loci = n_loci
  )
  loo_n <- (sum(n_sum) - n_sum) / (sum(n_loci) - n_loci)
  expect_equal(round(sum(n_sum) / sum(n_loci)), 19)
  expect_equal(round(loo_n), rep(c(18, 19), each = 5))

  jack <- .jackknife_tajimas_d_blocks(term_blocks)
  expect_equal(jack$estimate, tajimas_d_from_term_sums(colSums(term_blocks)))
  # the sample sizes of the blocks only matter through their total: the result
  # is the same when the total is spread evenly across blocks (so that
  # leaving out a block never changes the rounding)
  even_blocks <- term_blocks
  even_blocks[, "n"] <- sum(n_sum) / length(n_sum)
  expect_equal(.jackknife_tajimas_d_blocks(even_blocks), jack)
})

test_that("jackknife_pop_tajimas_d returns the right format", {
  test_gt <- load_example_gt("gen_tbl_gen_dist")
  jack_matrix <- jackknife_pop_tajimas_d(test_gt, block_size = 2)
  jack_tidy <- jackknife_pop_tajimas_d(test_gt, type = "tidy", block_size = 2)
  jack_list <- jackknife_pop_tajimas_d(test_gt, type = "list", block_size = 2)
  expect_equal(names(jack_matrix), c("quantity", "pop1", "pop2", "pop3"))
  expect_equal(
    names(jack_tidy),
    c("group", "estimate", "jack_estimate", "se", "z", "n_blocks")
  )
  expect_equal(names(jack_list), c("pop1", "pop2", "pop3"))
  expect_equal(jack_tidy$se, unname(unlist(jack_matrix[3, -1])))
  expect_equal(jack_list$pop2$se, jack_tidy$se[2])
  # ungrouped data returns a single row
  jack_ungrouped <- jackknife_pop_tajimas_d(dplyr::ungroup(test_gt),
    block_size = 2
  )
  expect_equal(nrow(jack_ungrouped), 1)
  expect_equal(
    jack_ungrouped$estimate,
    pop_tajimas_d(dplyr::ungroup(test_gt))
  )
})

test_that("jackknife_pop_tajimas_d works with several grouping variables", {
  test_gt <- load_example_gt("gen_tbl_gen_dist")
  test_gt$sex <- rep(c("F", "M"), length.out = nrow(test_gt))
  test_gt$pop_sex <- paste(test_gt$population, test_gt$sex, sep = ".")
  # the same results as grouping by a single variable combining them
  expect_equal(
    jackknife_pop_tajimas_d(
      dplyr::group_by(test_gt, population, sex),
      block_size = 2
    ),
    jackknife_pop_tajimas_d(dplyr::group_by(test_gt, pop_sex), block_size = 2)
  )
})

test_that("jackknife_pop_tajimas_d counts pseudohaploid allele copies", {
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

  jack_taj <- jackknife_pop_tajimas_d(
    pseudohaploid_gt,
    type = "list",
    block_size = 1
  )
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
  expect_equal(
    unname(vapply(jack_taj, function(x) x$estimate, numeric(1))),
    unname(unlist(expected))
  )
})

test_that("jackknife_pop_tajimas_d matches pop_tajimas_d with missing pi", {
  test_gt <- load_example_gt("gen_tbl_gen_dist")
  # make the 2nd locus missing for all individuals in pop1
  test_genotypes <- show_genotypes(test_gt)
  test_genotypes[dplyr::group_indices(test_gt) == 1, 2] <- NA
  test_gt_miss <- gen_tibble(
    x = test_genotypes,
    loci = show_loci(test_gt) %>% dplyr::select(-"big_index"),
    indiv_meta = test_gt %>%
      dplyr::ungroup() %>%
      dplyr::select(dplyr::all_of(c("id", "population"))),
    quiet = TRUE
  ) %>%
    dplyr::group_by(population)
  jack_tidy <- jackknife_pop_tajimas_d(test_gt_miss,
    type = "tidy", block_size = 2
  )
  expect_equal(jack_tidy$estimate, unlist(pop_tajimas_d(test_gt_miss)))
})
