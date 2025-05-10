test_that("gen_tibble_vcf filters non-biallelic SNPs", {
  # get path for non-biallelic SNPs
  vcf_path <- test_path("testdata/pop_a_non-biallelic.vcf")
  # read in the vcf
  vcf_cpp_gt <- gen_tibble_vcf(
    x = vcf_path,
    parser = "cpp",
    backingfile = tempfile(),
    quiet = TRUE
  )
  
})