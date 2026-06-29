skip_if_not_installed("vcfR")

# limit number of threads for tests
data.table::setDTthreads(2)
if (rlang::is_installed("RhpcBLASctl")) {
  RhpcBLASctl::blas_set_num_threads(2)
  RhpcBLASctl::omp_set_num_threads(2)
}

test_that("count_vcf_individuals with gzfile", {
  vcf_path <- system.file(
    "extdata/ploidy/ploidy_test.vcf.gz",
    package = "tidypopgen"
  )
  n_individuals <- count_vcf_individuals(vcf_path)
  vcfr_obj <- vcfR::read.vcfR(vcf_path, verbose = FALSE)
  expect_equal(n_individuals, ncol(vcfR::extract.gt(vcfr_obj)))
})
