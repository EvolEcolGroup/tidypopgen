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
  test_fbm <- tidypopgen:::.gt_get_bigsnp(test_gt)$genotypes
  test_as <- snp_allele_sharing(test_fbm)
  # use hierfstat to compute it
  dos_hier_match <- hierfstat::matching(test_genotypes)
  expect_true(all.equal(test_as, dos_hier_match))

  # check that we get the same result if we split the operation into two blocks
  test_as_2blocks <- snp_allele_sharing(test_fbm, block.size = 3)
  expect_identical(test_as_2blocks[], test_as[])

  # now estimate it with gen_tibble
  test_as_gt <- pairwise_allele_sharing(test_gt, as_matrix = TRUE)
  expect_true(all.equal(test_as, test_as_gt, check.attributes = FALSE))
})
