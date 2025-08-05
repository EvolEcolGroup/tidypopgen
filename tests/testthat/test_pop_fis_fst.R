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


test_that("pop_fst and pop_fist WG17 compute correctly", {
  # expect error if the tibble is not grouped
  expect_error(
    test_gt %>% pop_fis(method = "WG17"),
    ".x should be a grouped gen_tibble"
  )

  expect_error(
    test_gt %>% pop_fst(),
    ".x should be a grouped gen_tibble"
  )

  # group the tibble
  test_gt <- test_gt %>% dplyr::group_by(population)

  # compare results against raw hierfstat code
  fis_by_pop <- test_gt %>% pop_fis(include_global = TRUE, method = "WG17")
  fis_by_pop_hier <- hierfstat::fis.dosage(
    test_genotypes,
    pop = test_indiv_meta$population
  )
  expect_true(all.equal(fis_by_pop, fis_by_pop_hier, check.attributes = FALSE))
  # now check that we don't get the global
  fis_by_pop_sub <- test_gt %>% pop_fis(method = "WG17")
  expect_true(all.equal(fis_by_pop[-length(fis_by_pop)], fis_by_pop_sub))

  # compare results against raw hierfstat code
  fst_by_pop <- test_gt %>% pop_fst(include_global = TRUE)
  fst_by_pop_hier <- hierfstat::fst.dosage(
    test_genotypes,
    pop = test_indiv_meta$population
  )
  expect_true(all.equal(fst_by_pop, fst_by_pop_hier, check.attributes = FALSE))
  # now check that we don't get the global
  fst_by_pop_sub <- test_gt %>% pop_fst()
  expect_true(all.equal(fst_by_pop[-length(fst_by_pop)], fst_by_pop_sub))
})
