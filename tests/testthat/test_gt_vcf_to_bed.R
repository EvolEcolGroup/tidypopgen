testthat::test_that("as_gen_tibble works on genlight",{
  vcf_path <- system.file("/extdata/anolis/punctatus_t70_s10_n46_filtered.recode.vcf.gz",
                          package = "tidypopgen")
  bed_path <- gt_vcf_to_bed(vcf_path, bed_path = tempfile("anolis_"))
  test_gt <- gen_tibble(bed_path,quiet = TRUE)

  vcf_test <- vcfR::read.vcfR(file = vcf_path, verbose = FALSE)
  vcf_bi <- vcf_test[vcfR::is.biallelic(vcf_test)]
  expect_true(nrow(show_loci(test_gt))==nrow(vcf_bi@fix))
  expect_true(nrow(test_gt)==(ncol(vcf_bi@gt)-1))
})

