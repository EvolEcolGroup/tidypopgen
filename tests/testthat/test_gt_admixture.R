test_that("gt_admixture", {
  bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
  res_admix <- gt_admixture(bed_file, k = 3)

})
