skip_if_not_installed("hierfstat")

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

  # compare to hierfstat (we don't give L, so it is just the sum of loci pi)
  taj_hierfstat <- hierfstat::TajimaD.dosage(test_genotypes)
  expect_equal(taj_gt, taj_hierfstat)


  # repeat the tests for a subset of the data
  # remove the 2nd individual and the 3rd and 5th snp
  test_genotypes_subset1 <- test_genotypes[-2, c(-3, -5)]
  test_gt_subset1 <- test_gt %>%
    filter(id != "b") %>%
    select_loci(c(-3, -5))
  taj_gt <- pop_tajimas_d(test_gt_subset1$genotypes)
  taj_hierfstat <- hierfstat::TajimaD.dosage(test_genotypes_subset1)
  expect_equal(taj_gt, taj_hierfstat)
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
  ploidy_by_group <- pseudohaploid_gt %>%
    dplyr::summarise(n = sum(indiv_ploidy(genotypes))) %>%
    dplyr::pull(.data$n)

  expected <- Map(tajimas_d_from_pi_vec, pi_by_group, ploidy_by_group)
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
  expect_equal(taj_gt, tajimas_d_from_pi_vec(loci_pi(test_gt), 6 * 4))
  # grouped
  test_gt <- test_gt %>% group_by(population)
  taj_map <- test_gt %>% group_map(.f = ~ pop_tajimas_d(.x))
  taj_grp <- test_gt %>% pop_tajimas_d(n_cores = 2)
  expect_equal(taj_map, taj_grp)
  pi_by_group <- loci_pi(test_gt, type = "list")
  expect_equal(
    unlist(taj_grp),
    unname(unlist(Map(tajimas_d_from_pi_vec, pi_by_group, c(12, 12))))
  )
})
