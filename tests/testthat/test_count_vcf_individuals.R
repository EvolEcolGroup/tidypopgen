test_that("count_vcf_individuals with gzfile",{
  vcf_path <- system.file("/extdata/anolis/punctatus_t70_s10_n46_filtered.recode.vcf.gz",
                          package = "tidypopgen")
  n_individuals <- count_vcf_individuals(vcf_path)
  vcfr_obj <- vcfR::read.vcfR(vcf_path, verbose = FALSE)
  expect_true(n_individuals == dim(vcfr_obj)[2])
})