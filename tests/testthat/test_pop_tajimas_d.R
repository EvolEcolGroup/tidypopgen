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
