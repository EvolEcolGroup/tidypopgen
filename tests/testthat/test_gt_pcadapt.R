test_that("gt_pcadapt works on gt_pca object",{
  bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
  missing_gt <- gen_tibble(bed_file,  backingfile = tempfile("missing_"),quiet=TRUE)
  missing_gt <- gt_impute_simple(missing_gt, method = "mode")
  missing_pca <- missing_gt %>% gt_pca_partialSVD()
  missing_pcadapt <- gt_pcadapt(missing_gt, missing_pca, k = 3)
  expect_true((inherits(missing_pcadapt, "gt_pcadapt")))
})


