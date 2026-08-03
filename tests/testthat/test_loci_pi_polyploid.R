# tests for loci_pi() on polyploid and mixed-ploidy gen_tibbles.
# hierfstat::pi.dosage() (used to validate the diploid case in
# test_loci_pi.R) is diploid-only, so here we validate against an
# independent reference implementation computed directly from the raw
# dosage matrix.

# reference implementation: c_0 * c_1 / (M * (M-1) / 2), where M is the
# total number of valid allele copies (sum of ploidy over non-missing
# individuals), computed per locus (optionally within each population).
ref_loci_pi <- function(genotypes, ploidy) {
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

ref_grouped_loci_pi <- function(genotypes, ploidy, population) {
  pop_levels <- unique(population)
  out <- sapply(pop_levels, function(pop) { # nolint
    rows <- population == pop
    ref_loci_pi(genotypes[rows, , drop = FALSE], ploidy[rows])
  })
  colnames(out) <- pop_levels
  out
}

test_that("loci_pi computes correctly for uniform tetraploid data", {
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
  expect_true(all(indiv_ploidy(test_gt) == 4))

  expected <- ref_loci_pi(test_genotypes, test_ploidy)
  expect_equal(loci_pi(test_gt), unname(expected))
  expect_equal(loci_pi(test_gt$genotypes), unname(expected))
})

test_that("loci_pi computes correctly for mixed-ploidy data", {
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
  expect_identical(indiv_ploidy(test_gt), test_ploidy)
  expect_true(show_ploidy(test_gt) == 0)

  expected <- ref_loci_pi(test_genotypes, test_ploidy)
  expect_equal(loci_pi(test_gt), unname(expected))
})

test_that("loci_pi computes correctly grouped by uniform-ploidy populations", { # nolint
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

  expected <- ref_grouped_loci_pi(
    test_genotypes,
    test_ploidy,
    test_indiv_meta$population
  )

  test_pi <- test_gt %>% loci_pi(type = "matrix")
  expect_equal(unname(test_pi), unname(expected))
  expect_identical(colnames(test_pi), c("popA", "popB"))
})

test_that("loci_pi handles a population mixing diploid and tetraploid individuals", { # nolint
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

  expected <- ref_grouped_loci_pi(
    test_genotypes,
    test_ploidy,
    test_indiv_meta$population
  )
  test_pi <- test_gt %>% loci_pi(type = "matrix")
  expect_equal(unname(test_pi), unname(expected))
})

test_that("loci_pi is NA when a locus has 0 or 1 valid allele copies in a group", { # nolint
  test_indiv_meta <- data.frame(
    id = letters[1:3],
    population = c("popA", "popA", "popB")
  )
  test_genotypes <- rbind(
    c(1, NA),
    c(2, NA),
    c(NA, 1)
  )
  test_ploidy <- c(4, 4, 1)
  test_loci <- data.frame(
    name = paste0("rs", 1:2),
    chromosome = c(1, 1),
    position = c(3, 5),
    genetic_dist = as.double(rep(0, 2)),
    allele_ref = c("A", "T"),
    allele_alt = c("T", "C")
  )

  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    ploidy = test_ploidy,
    quiet = TRUE
  ) %>% dplyr::group_by(population)

  test_pi <- test_gt %>% loci_pi(type = "matrix")
  # popB, locus 2: 0 valid individuals -> NA
  expect_true(is.na(test_pi["rs2", "popB"]))
  # popB, locus 1: a single haploid individual -> M = 1 -> NA
  expect_true(is.na(test_pi["rs1", "popB"]))
})

test_that("the polyploid kernel reduces to the diploid kernel at ploidy 2", {
  test_genotypes <- rbind(
    c(1, 1, 0, 1, 1, 2),
    c(2, 1, 0, NA, 0, NA),
    c(2, 2, 0, 0, 1, NA),
    c(1, 0, 0, 1, 2, 0)
  )
  ploidy_2 <- rep(2, nrow(test_genotypes))

  test_indiv_meta <- data.frame(
    id = letters[1:4],
    population = c("pop1", "pop1", "pop2", "pop2")
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
    BM = X, rowInd = rows, colInd = ind, ploidy = ploidy_2, ncores = 1
  )
  expect_equal(dip_pi, poly_pi)

  grouped <- test_gt %>% dplyr::group_by(population)
  dip_grouped <- tidypopgen:::gt_grouped_pi_diploid(
    BM = X, rowInd = rows, colInd = ind,
    groupIds = dplyr::group_indices(grouped) - 1,
    ngroups = 2, ncores = 1
  )$pi
  poly_grouped <- tidypopgen:::gt_grouped_pi_polyploid(
    BM = X, rowInd = rows, colInd = ind,
    groupIds = dplyr::group_indices(grouped) - 1,
    ngroups = 2, ploidy = ploidy_2, ncores = 1
  )$pi
  expect_equal(dip_grouped, poly_grouped)
})
