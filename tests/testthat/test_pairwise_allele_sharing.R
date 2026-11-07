skip_if_not_installed("hierfstat")

test_genotypes <- rbind(
  c(1, 1, 0, 1, 1, 0),
  c(2, 1, 0, NA, 0, 0),
  c(2, NA, 0, 0, 1, 1)
)
test_indiv_meta <- data.frame(
  id = c("a", "b", "c"),
  population = c("pop1", "pop1", "pop2")
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


test_that("snp_allele_sharing and pairwise_allele_sharing computes allele sharing correctly", { # nolint
  test_fbm <- tidypopgen:::.gt_get_fbm(test_gt)
  test_as <- snp_allele_sharing(test_fbm)
  # use hierfstat to compute it
  dos_hier_match <- hierfstat::matching(test_genotypes)
  expect_true(all.equal(test_as, dos_hier_match))

  # check that we get the same result if we split the operation into two blocks
  test_as_2blocks <- snp_allele_sharing(test_fbm, block.size = 3)
  expect_identical(test_as_2blocks[], test_as[])

  # now estimate it with gen_tibble
  test_as_gt <- pairwise_allele_sharing(test_gt, as_matrix = TRUE)
  # remove the special class to allow comparison with the reference
  class(test_as_gt) <- setdiff(class(test_as_gt), "pairwise_matrix")
  expect_true(all.equal(test_as, test_as_gt,
    check.attributes = FALSE
  ))
})

test_that("snp_allele_sharing is correct when blocks have different sizes", {
  test_fbm <- tidypopgen:::.gt_get_fbm(test_gt)
  # 5 loci in blocks of at most 2 loci gives blocks of different sizes
  test_as_uneven <- snp_allele_sharing(test_fbm,
    ind.col = 1:5,
    block.size = 2
  )
  dos_hier_match <- hierfstat::matching(test_genotypes[, 1:5])
  expect_equal(test_as_uneven, dos_hier_match, ignore_attr = TRUE)
  # and the same as computing it in a single block
  expect_equal(
    test_as_uneven,
    snp_allele_sharing(test_fbm, ind.col = 1:5)
  )
})

test_that("snp_allele_sharing is correct across bit-packing boundaries", {
  # more than 64 individuals and loci, with missing data, so that genotypes
  # are packed into several words and blocks of individuals
  set.seed(123)
  geno_vals <- sample(c(0:2, NA), 70 * 150,
    replace = TRUE,
    prob = c(0.4, 0.3, 0.25, 0.05)
  )
  geno_mat <- matrix(geno_vals, nrow = 70)
  big_gt <- gen_tibble(
    x = geno_mat,
    loci = data.frame(
      name = paste0("rs", 1:150), chromosome = 1L, position = 1:150,
      genetic_dist = 0, allele_ref = "A", allele_alt = "T"
    ),
    indiv_meta = data.frame(id = paste0("i", 1:70), population = "p1"),
    quiet = TRUE
  )
  test_fbm <- tidypopgen:::.gt_get_fbm(big_gt)
  # reordered subsets of individuals and loci, and multiple cores
  rows <- c(70:3, 1)
  cols <- c(150:80, 2:60)
  test_as <- snp_allele_sharing(test_fbm,
    ind.row = rows, ind.col = cols,
    block.size = 70, n_cores = 2
  )
  # reference: 1/2 * (1 + tcrossprod(dos - 1) / tcrossprod(valid)), with
  # missing dosages set to 1
  g <- geno_mat[rows, cols]
  valid <- !is.na(g) + 0
  dos_1 <- g - 1
  dos_1[is.na(dos_1)] <- 0
  as_ref <- 0.5 * (1 + tcrossprod(dos_1) / tcrossprod(valid))
  expect_equal(test_as, as_ref, tolerance = 1e-12)
})

test_that("pairwise_allele_sharing uses imputed genotypes only if set to", {
  # test_gt has missing genotypes
  as_na <- pairwise_allele_sharing(test_gt, as_matrix = TRUE)
  gt_imputed <- gt_impute_simple(test_gt, method = "mode")
  # imputed genotypes are missing unless the gen_tibble is set to use them
  gt_set_imputed(gt_imputed, set = FALSE)
  expect_identical(
    pairwise_allele_sharing(gt_imputed, as_matrix = TRUE),
    as_na
  )
  # when used, they are the same as genotypes that were never missing
  gt_set_imputed(gt_imputed, set = TRUE)
  gt_filled <- gen_tibble(
    x = show_genotypes(gt_imputed), loci = test_loci,
    indiv_meta = test_indiv_meta, quiet = TRUE
  )
  expect_identical(
    pairwise_allele_sharing(gt_imputed, as_matrix = TRUE),
    pairwise_allele_sharing(gt_filled, as_matrix = TRUE)
  )
})
