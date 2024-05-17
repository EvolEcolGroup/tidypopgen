test_that("import a vcf with multiple ploidy",{
  vcf_path <- system.file("/extdata/ploidy_test.vcf.gz",
                          package = "tidypopgen")
  test_gt <- gen_tibble(vcf_path, backingfile = tempfile(), quiet = TRUE)
  # it's mixed ploidy
  expect_true(show_ploidy(test_gt)==0)
  # individuals are either 2 or 4
  expect_true(all(indiv_ploidy(test_gt) %in% c(2,4)))
})
