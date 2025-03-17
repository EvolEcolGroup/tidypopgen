bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
missing_gt <- gen_tibble(
  bed_file,
  backingfile = tempfile("missing_"),
  quiet = TRUE
)
missing_gt <- gt_impute_simple(missing_gt, method = "mode")
missing_pca <- missing_gt %>% gt_pca_partialSVD()

test_that("gt_pcadapt works on gt_pca object", {
  missing_pcadapt <- gt_pcadapt(missing_gt, missing_pca, k = 3)
  expect_true((inherits(missing_pcadapt, "gt_pcadapt")))
})

test_that("n_cores can be set", {
  expect_true(getOption("bigstatsr.check.parallel.blas"))
  one_core <- gt_pcadapt(missing_gt, missing_pca, k = 3, n_cores = 1)
  expect_true(getOption("bigstatsr.check.parallel.blas"))
  two_core <- gt_pcadapt(missing_gt, missing_pca, k = 3, n_cores = 2)
  expect_true(getOption("bigstatsr.check.parallel.blas"))
})
