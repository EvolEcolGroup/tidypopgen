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
  # @TODO why does the above have names with one parser and NOT the other?!?
  expect_true(all.equal(show_ploidy(test_gt), show_ploidy(test_cpp_gt)))
  expect_true(all.equal(show_loci(test_gt), show_loci(test_cpp_gt)))
})
