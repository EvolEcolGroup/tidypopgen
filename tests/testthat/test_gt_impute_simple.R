test_that("impute and use the imputation",{
  bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
  missing_gt <- gen_tibble(bed_file,  backingfile = tempfile("missing_"), quiet = TRUE)
  # we get errors because of missing values
  expect_error(missing_gt %>% gt_pca_partialSVD(),
               "You can't have missing values")
  expect_false(gt_has_imputed(missing_gt))
  expect_error(gt_uses_imputed(missing_gt),
               "this dataset does not have any imputated")
  # now impute
  missing_gt <- gt_impute_simple(missing_gt)
  # we have imputed
  expect_true(gt_has_imputed(missing_gt))
  # but don't use it by default
  expect_false(gt_uses_imputed(missing_gt))
  # now we return a pca successfully
  expect_true(inherits(missing_gt %>% gt_pca_partialSVD(),"gt_pca"))
  # simple error message
  expect_error(gt_set_imputed(missing_gt),
               "set should be either")
})
