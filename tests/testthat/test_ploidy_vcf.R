skip_if_not_installed("vcfR")

test_that("import a vcf with multiple ploidy", {
  vcf_path <- system.file(
    "/extdata/ploidy/ploidy_test.vcf.gz",
    package = "tidypopgen"
  )
  test_gt <- gen_tibble(vcf_path,
    backingfile = tempfile(), quiet = TRUE,
    parser = "vcfR"
  )
  # it's mixed ploidy
  expect_true(show_ploidy(test_gt) == 0)
  # individuals are either 2 or 4
  expect_true(all(indiv_ploidy(test_gt) %in% c(2, 4)))
  # now try the cpp parser
  test_cpp_gt <- gen_tibble(
    vcf_path,
    backingfile = tempfile(),
    quiet = TRUE,
    parser = "cpp"
  )
  expect_true(all.equal(show_genotypes(test_gt), show_genotypes(test_cpp_gt)))
  expect_true(all.equal(indiv_ploidy(test_gt), indiv_ploidy(test_cpp_gt)))
  expect_true(all.equal(show_ploidy(test_gt), show_ploidy(test_cpp_gt)))
  expect_true(all.equal(show_loci(test_gt), show_loci(test_cpp_gt)))
})

test_that("indiv_het_obs and pop_het_obs run end-to-end on real mixed-ploidy data", { # nolint
  vcf_path <- system.file(
    "/extdata/ploidy/ploidy_test.vcf.gz",
    package = "tidypopgen"
  )
  test_gt <- gen_tibble(vcf_path,
    backingfile = tempfile(), quiet = TRUE,
    parser = "vcfR"
  )
  # sanity check on the fixture itself: a genuine mix of diploid/tetraploid
  expect_true(all(indiv_ploidy(test_gt) %in% c(2, 4)))
  expect_true(length(unique(indiv_ploidy(test_gt))) == 2)

  het_obs <- indiv_het_obs(test_gt)
  expect_true(all(het_obs >= 0 & het_obs <= 1, na.rm = TRUE))

  # group individuals by their id prefix (encodes both site and ploidy,
  # e.g. "ta" = tetraploid, "da" = diploid) to get within-ploidy population
  # groups, and one group ("BAL") that also works when ploidy is uniform
  test_gt$population <- sub("_[0-9]+.*$", "", test_gt$id)
  test_gt_grp <- test_gt %>% dplyr::group_by(population)

  pop_ho <- test_gt_grp %>% pop_het_obs(by_locus = TRUE)
  expect_true(all(pop_ho >= 0 & pop_ho <= 1, na.rm = TRUE))
  expect_identical(colnames(pop_ho), sort(unique(test_gt$population)))

  pop_ho_global <- test_gt_grp %>%
    pop_het_obs(by_locus = TRUE, include_global = TRUE)
  expect_equal(
    unname(pop_ho_global[, "global"]),
    rowMeans(pop_ho, na.rm = TRUE)
  )
})
