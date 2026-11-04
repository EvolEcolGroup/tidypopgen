# limit number of threads for tests
data.table::setDTthreads(2)
if (rlang::is_installed("RhpcBLASctl")) {
  RhpcBLASctl::blas_set_num_threads(2)
  RhpcBLASctl::omp_set_num_threads(2)
}

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


test_that("indiv_inbreeding compute correctly", {
  # expect error if the tibble is not grouped
  inbreed_hier <- diag(hierfstat::beta.dosage(test_genotypes))
  inbreed_gt <- test_gt %>% indiv_inbreeding(method = "WG17")
  expect_true(all.equal(inbreed_hier, inbreed_gt, check.attributes = FALSE))

  # group the tibble
  test_gt <- test_gt %>% dplyr::group_by(population)
  # is this right???
  inbreed_group_gt <- test_gt %>% indiv_inbreeding(method = "WG17")
  # subset gentibble to a single population, pop 1
  pop1_gt <- test_gt %>% dplyr::filter(population == "pop1")
  inbreed_pop1_hier <- diag(hierfstat::beta.dosage(show_genotypes(pop1_gt)))
  expect_true(all.equal(
    inbreed_pop1_hier,
    inbreed_group_gt[[1]],
    check.attributes = FALSE
  ))
  # subset gentibble to a single population, pop 2
  pop2_gt <- test_gt %>% dplyr::filter(population == "pop2")
  inbreed_pop2_hier <- diag(hierfstat::beta.dosage(show_genotypes(pop2_gt)))
  expect_true(all.equal(
    inbreed_pop2_hier,
    inbreed_group_gt[[2]],
    check.attributes = FALSE
  ))
})

test_that("indiv_inbreeding and allele sharing work with a single individual", {
  test_gt <- load_example_gt("gen_tbl")
  single_gt <- test_gt[1, ]
  # a 1 x 1 matrix, rather than a vector
  ashare <- pairwise_allele_sharing(single_gt, as_matrix = TRUE)
  expect_equal(dim(ashare), c(1, 1))
  expect_equal(dimnames(ashare), list(single_gt$id, single_gt$id))
  # the same value as in the matrix of all individuals
  expect_equal(
    unclass(ashare)[1, 1],
    unclass(pairwise_allele_sharing(test_gt, as_matrix = TRUE))[1, 1]
  )
  expect_equal(dim(pairwise_grm(single_gt)), c(1, 1))
  # inbreeding is relative to the other individuals, so it is NaN for a
  # single individual
  expect_true(is.nan(indiv_inbreeding(single_gt)))
  # a group with a single individual
  test_gt$population[1] <- "pop_single"
  grouped_gt <- test_gt %>% group_by(population)
  inb_grp <- indiv_inbreeding(grouped_gt)
  expect_length(inb_grp, 4)
  single_grp <- which(group_keys(grouped_gt)$population == "pop_single")
  expect_true(is.nan(inb_grp[[single_grp]]))
})
