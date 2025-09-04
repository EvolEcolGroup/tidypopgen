skip_if_not_installed("hierfstat")

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

test_that("pop_* basic stats functions work correctly", {
  test_gt <- test_gt %>% dplyr::group_by(population)
  test_hier <- gt_as_hierfstat(test_gt)
  basic_hier <- hierfstat::basic.stats(test_hier)

  # observed heterozygosity by locus
  test_het_obs <- test_gt %>% pop_het_obs(by_locus = TRUE)
  expect_true(all.equal(
    basic_hier$Ho,
    round(test_het_obs, 4),
    check.attributes = FALSE
  ))
  # overall (mean of by locus values)
  test_het_obs <- test_gt %>% pop_het_obs(by_locus = FALSE)
  expect_true(all.equal(
    colMeans(basic_hier$Ho, na.rm = TRUE),
    test_het_obs,
    check.attributes = FALSE
  ))
  # check the overall value
  test_het_obs <- test_gt %>%
    pop_het_obs(by_locus = TRUE, include_global = TRUE)
  expect_true(all.equal(
    round(test_het_obs[, 4], 4),
    basic_hier$perloc$Ho,
    check.attributes = FALSE
  ))
  test_het_obs <- test_gt %>%
    pop_het_obs(by_locus = FALSE, include_global = TRUE)
  expect_true(all.equal(
    round(test_het_obs[4], 4),
    basic_hier$overall["Ho"],
    check.attributes = FALSE
  ))

  # test expected heterozygosity by locus
  test_het_exp <- test_gt %>% pop_het_exp(by_locus = TRUE)
  expect_true(all.equal(
    basic_hier$Hs,
    round(test_het_exp, 4),
    check.attributes = FALSE
  ))
  # overall (mean of by locus values)
  test_het_exp <- test_gt %>% pop_het_exp(by_locus = FALSE)
  expect_true(all.equal(
    colMeans(basic_hier$Hs, na.rm = TRUE),
    test_het_exp,
    check.attributes = FALSE
  ))
  # check the overall value
  test_het_exp <- test_gt %>%
    pop_het_exp(by_locus = TRUE, include_global = TRUE)
  expect_true(all.equal(
    round(test_het_exp[, 4], 4),
    basic_hier$perloc$Hs,
    check.attributes = FALSE
  ))
  test_het_exp <- test_gt %>%
    pop_het_exp(by_locus = FALSE, include_global = TRUE)
  expect_true(all.equal(
    round(test_het_exp[4], 4),
    basic_hier$overall["Hs"],
    check.attributes = FALSE
  ))

  # test Fis by locus
  test_fis <- test_gt %>% pop_fis(by_locus = TRUE, method = "Nei87")
  expect_true(all.equal(
    basic_hier$Fis,
    round(test_fis, 4),
    check.attributes = FALSE
  ))
  # overall (mean of by locus values)
  test_fis <- test_gt %>% pop_fis(by_locus = FALSE)
  expect_true(all.equal(
    round(colMeans(basic_hier$Fis, na.rm = TRUE), 4),
    round(test_fis, 4),
    check.attributes = FALSE
  ))
  # check the overall value
  test_fis <- test_gt %>% pop_fis(by_locus = TRUE, include_global = TRUE)
  expect_true(all.equal(
    round(test_fis[, 4], 4),
    basic_hier$perloc$Fis,
    check.attributes = FALSE
  ))
  test_fis <- test_gt %>% pop_fis(by_locus = FALSE, include_global = TRUE)
  expect_true(all.equal(
    round(test_fis[4], 4),
    basic_hier$overall["Fis"],
    check.attributes = FALSE
  ))

  # test global stats by locus
  test_global_stats <- test_gt %>% pop_global_stats(by_locus = TRUE)
  expect_true(all.equal(
    basic_hier$perloc,
    round(test_global_stats, 4),
    check.attributes = FALSE
  ))
  # test global stats overall
  test_global_stats <- test_gt %>% pop_global_stats(by_locus = FALSE)
  expect_true(all.equal(
    basic_hier$overall,
    round(test_global_stats, 4),
    check.attributes = FALSE
  ))

  # test Fis errors for incorrect parameters for Nei87
  expect_error(
    test_gt %>%
      pop_fis(
        by_locus = TRUE,
        method = "Nei87",
        allele_sharing_mat = matrix(1, nrow = 7, ncol = 7)
      ),
    "allele_sharing_mat not relevant for Nei87"
  )
  # test Fis errors for incorrect parameters for WG17
  expect_error(
    test_gt %>%
      pop_fis(
        by_locus = TRUE,
        method = "WG17",
        allele_sharing_mat = matrix(1, nrow = 7, ncol = 7)
      ),
    "by_locus not implemented for WG17"
  )
  # test alias for pop_het_exp
  expect_identical(
    test_het_exp,
    test_gt %>% pop_gene_div(by_locus = FALSE, include_global = TRUE)
  )
})

test_that("pop_* basic stats work on a single population", {
  test_het_obs <- test_gt %>% pop_het_obs(by_locus = TRUE)
  # just one column, as we only have one population
  expect_true(ncol(test_het_obs) == 1)
  test_fis <- test_gt %>% pop_fis(by_locus = TRUE)
  # just one column, as we only have one population
  expect_true(ncol(test_fis) == 1)
  test_het_exp <- test_gt %>% pop_het_exp(by_locus = TRUE)
  # just one column, as we only have one population
  expect_true(ncol(test_het_exp) == 1)
  test_global_stats <- test_gt %>% pop_global_stats(by_locus = TRUE)
  expect_true(all(is.nan(test_global_stats$Fstp)))
  # Fstp is not defined for a single population
})

test_that("n_cores can be set", {
  # pop_het_obs
  one_core <- test_gt %>% pop_het_obs(by_locus = TRUE, n_cores = 1)
  two_core <- test_gt %>% pop_het_obs(by_locus = TRUE, n_cores = 2)
  expect_equal(one_core, two_core)
  # pop_het_exp
  one_core <- test_gt %>% pop_het_exp(by_locus = TRUE, n_cores = 1)
  two_core <- test_gt %>% pop_het_exp(by_locus = TRUE, n_cores = 2)
  expect_equal(one_core, two_core)
  # pop_global_stats
  one_core <- test_gt %>% pop_global_stats(by_locus = TRUE, n_cores = 1)
  two_core <- test_gt %>% pop_global_stats(by_locus = TRUE, n_cores = 2)
  expect_equal(one_core, two_core)
})
