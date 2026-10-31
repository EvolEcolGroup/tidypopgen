# Tests of how gt_admixture() chooses its copy of ADMIXTURE, which do not need
# ADMIXTURE to be installed

test_that("without conda, conda_env = 'auto' uses a copy on the PATH", {
  skip_if_not_installed("reticulate")
  bed_file <- system.file("extdata", "lobster", "lobster.bed",
    package = "tidypopgen"
  )
  # conda_list() gives an error when conda itself is not installed
  local_mocked_bindings(
    conda_list = function(...) {
      stop("Unable to find conda binary. Is Anaconda installed?")
    },
    .package = "reticulate"
  )
  res <- tryCatch(
    gt_admixture(bed_file, k = 2, seed = 1),
    error = function(e) conditionMessage(e)
  )
  # the copy of ADMIXTURE on the PATH is used, or its absence is reported,
  # rather than the error of conda_list()
  if (nzchar(Sys.which("admixture"))) {
    expect_s3_class(res, "gt_admix")
  } else {
    expect_identical(res, "ADMIXTURE is not installed in this path")
  }
})
