test_that("gt_admixture", {
  # skip if admixture is not installed (TODO we should also check if we have it through reticulate)
  skip_if ((system2("which", args = "admixture") != 0))

  bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
  res_admix <- gt_admixture(bed_file, k = 3)

})
