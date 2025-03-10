test_that("count_vcf_variants with gzfile", {
  vcf_path <- system.file("/extdata/ploidy/ploidy_test.vcf.gz",
    package = "tidypopgen"
  )
  n_variants <- count_vcf_variants(vcf_path)
  vcfr_obj <- vcfR::read.vcfR(vcf_path, verbose = FALSE)
  expect_true(n_variants == dim(vcfr_obj)[1])
})
