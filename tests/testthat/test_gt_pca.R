test_that("fit_gt_pca_and_predict",{
  bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
  missing_gt <- gen_tibble(bed_file,  backingfile = tempfile("missing_"),quiet=TRUE)
  expect_error( missing_gt %>% gt_pca_partialSVD(),
                "You can't have missing")
  missing_gt <- gt_impute_simple(missing_gt, method = "mode")
  missing_pca <- missing_gt %>% gt_pca_partialSVD()
  expect_true(all.equal(predict(missing_pca),
                        predict(missing_pca, new_data = missing_gt),
                        check.attributes=FALSE))
})

