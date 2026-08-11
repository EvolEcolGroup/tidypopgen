# tests for pop_tajimas_d() and windows_pop_tajimas_d() on polyploid and
# mixed-ploidy gen_tibbles. hierfstat::TajimaD.dosage() (used to validate
# the diploid case in test_pop_tajimas_d.R) is diploid-only, so here we
# validate against an independent reference implementation computed
# directly from the raw dosage matrix.

# reference nucleotide diversity per locus, mirroring the reference in
# test_loci_pi_polyploid.R but reimplemented independently here
ref_pi <- function(genotypes, ploidy) {
  apply(genotypes, 2, function(col) {
    valid <- !is.na(col)
    if (!any(valid)) {
      return(NA_real_)
    }
    m_total <- sum(ploidy[valid])
    if (m_total <= 1) {
      return(NA_real_)
    }
    c_1 <- sum(col[valid])
    c_0 <- m_total - c_1
    c_0 * c_1 / (m_total * (m_total - 1) / 2)
  })
}

# reference Tajima's D, reimplemented independently from the package's
# internal tajimas_d_from_pi_vec helper
ref_tajimas_d <- function(pi, n) {
  seg <- sum(pi > 0 & pi < 1, na.rm = TRUE)
  k_hat <- sum(pi, na.rm = TRUE)
  a1 <- sum(1 / seq_len(n - 1))
  a2 <- sum(1 / seq_len(n - 1)^2)
  b1 <- (n + 1) / (3 * (n - 1))
  b2 <- 2 * (n^2 + n + 3) / (9 * n * (n - 1))
  c1 <- b1 - 1 / a1
  c2 <- b2 - (n + 2) / (a1 * n) + a2 / a1^2
  e1 <- c1 / a1
  e2 <- c2 / (a1^2 + a2)
  vd <- e1 * seg + e2 * seg * (seg - 1)
  (k_hat - seg / a1) / sqrt(vd)
}

test_that("pop_tajimas_d computes correctly for uniform tetraploid data", {
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
  )

  expected <- ref_tajimas_d(
    ref_pi(test_genotypes, test_ploidy),
    sum(test_ploidy)
  )
  expect_equal(pop_tajimas_d(test_gt$genotypes), expected)
  expect_equal(pop_tajimas_d(test_gt), expected)
})

test_that("pop_tajimas_d computes correctly for mixed-ploidy data", {
  test_indiv_meta <- data.frame(
    id = letters[1:4],
    population = rep("popA", 4)
  )
  # a,b diploid; c,d tetraploid
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
  )

  expected <- ref_tajimas_d(
    ref_pi(test_genotypes, test_ploidy),
    sum(test_ploidy)
  )
  expect_equal(pop_tajimas_d(test_gt), expected)
})

test_that("pop_tajimas_d computes correctly grouped by uniform-ploidy populations", { # nolint
  test_indiv_meta <- data.frame(
    id = letters[1:6],
    population = c("popA", "popA", "popA", "popB", "popB", "popB")
  )
  # popA tetraploid, popB diploid
  test_genotypes <- rbind(
    c(1, 4, 0, 2),
    c(4, 4, 0, 0),
    c(0, 2, 0, 4),
    c(1, 0, 2, 1),
    c(0, 2, 1, 0),
    c(2, 1, 0, NA)
  )
  test_ploidy <- c(4, 4, 4, 2, 2, 2)
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

  pop_a_rows <- 1:3
  pop_b_rows <- 4:6
  expected_a <- ref_tajimas_d(
    ref_pi(test_genotypes[pop_a_rows, ], test_ploidy[pop_a_rows]),
    sum(test_ploidy[pop_a_rows])
  )
  expected_b <- ref_tajimas_d(
    ref_pi(test_genotypes[pop_b_rows, ], test_ploidy[pop_b_rows]),
    sum(test_ploidy[pop_b_rows])
  )

  taj_grp <- test_gt %>% pop_tajimas_d(n_cores = 2)
  expect_equal(taj_grp[[1]], expected_a)
  expect_equal(taj_grp[[2]], expected_b)

  # matches group_map, mirroring the diploid grouped test
  taj_map <- test_gt %>% group_map(.f = ~ pop_tajimas_d(.x))
  expect_true(all.equal(taj_map, taj_grp))
})

test_that("pop_tajimas_d handles a population mixing diploid and tetraploid individuals", { # nolint
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

  expected <- ref_tajimas_d(
    ref_pi(test_genotypes, test_ploidy),
    sum(test_ploidy)
  )
  taj_grp <- test_gt %>% pop_tajimas_d()
  expect_equal(taj_grp[[1]], expected)
})

test_that("windows_pop_tajimas_d computes correctly for tetraploid data", {
  test_indiv_meta <- data.frame(
    id = letters[1:4],
    population = rep("popA", 4)
  )
  test_genotypes <- rbind(
    c(1, 4, 0, 2, 1, 3),
    c(4, 4, 0, 0, 2, 1),
    c(0, 2, 0, 4, 0, 4),
    c(2, 0, 4, NA, 3, 2)
  )
  test_ploidy <- rep(4, 4)
  test_loci <- data.frame(
    name = paste0("rs", 1:6),
    chromosome = rep(1, 6),
    position = c(3, 5, 23, 456, 500, 600),
    genetic_dist = as.double(rep(0, 6)),
    allele_ref = c("A", "T", "C", "G", "C", "T"),
    allele_alt = c("T", "C", "G", "A", "T", "A")
  )
  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    ploidy = test_ploidy,
    quiet = TRUE
  )

  window_taj <- windows_pop_tajimas_d(
    test_gt,
    window_size = 3,
    step_size = 3,
    size_unit = "snp",
    min_loci = 1
  )

  pi_by_locus <- ref_pi(test_genotypes, test_ploidy)
  expected_window1 <- ref_tajimas_d(pi_by_locus[1:3], sum(test_ploidy))
  expected_window2 <- ref_tajimas_d(pi_by_locus[4:6], sum(test_ploidy))

  expect_equal(window_taj$stat[1], expected_window1)
  expect_equal(window_taj$stat[2], expected_window2)
})

test_that("the polyploid pi kernel gives the same tajimas_d at ploidy 2 as the diploid kernel", { # nolint
  test_genotypes <- rbind(
    c(1, 1, 0, 1, 1, 2),
    c(2, 1, 0, NA, 0, NA),
    c(2, 2, 0, 0, 1, NA)
  )
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c"),
    population = c("pop1", "pop1", "pop2")
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

  X <- attr(test_gt$genotypes, "fbm") # nolint
  rows <- vctrs::vec_data(test_gt$genotypes)
  ind <- attr(test_gt$genotypes, "loci")$big_index

  dip_pi <- tidypopgen:::gt_pi_diploid(
    BM = X, rowInd = rows, colInd = ind, ncores = 1
  )
  poly_pi <- tidypopgen:::gt_pi_polyploid(
    BM = X, rowInd = rows, colInd = ind, ploidy = rep(2, 3), ncores = 1
  )
  expect_equal(dip_pi, poly_pi)

  n <- nrow(test_genotypes) * 2
  expect_equal(
    tidypopgen:::tajimas_d_from_pi_vec(dip_pi, n),
    tidypopgen:::tajimas_d_from_pi_vec(poly_pi, n)
  )
})
