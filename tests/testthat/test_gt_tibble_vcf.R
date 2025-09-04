skip_if_not_installed("vcfR")

test_that("gen_tibble_vcf filters non-biallelic SNPs", {
  # get path for non-biallelic SNPs
  vcf_path <- test_path("testdata/vcf", "pop_a_non_biallelic.vcf")
  # read in the vcf
  vcf_cpp_gt <- gen_tibble(
    x = vcf_path,
    parser = "cpp",
    backingfile = tempfile(),
    quiet = TRUE
  )
  pop_a_gt <-
    gen_tibble(
      system.file("extdata/pop_a.vcf", package = "tidypopgen"),
      backingfile = tempfile(),
      quiet = TRUE
    )
  # the vcf is missing locus rs3094315 and rs 307354
  expect_true(all(!c("rs3094315", "rs307354") %in% show_loci(vcf_cpp_gt)$name))
  # check that the genotypes are the same
  # remove them from the pop_a_gt tibble
  pop_a_gt_sub <- pop_a_gt %>% select_loci(!any_of(c("rs3094315", "rs307354")))
  # and check
  expect_true(all(show_genotypes(pop_a_gt_sub) == show_genotypes(vcf_cpp_gt)))

  # note that our vcfr does not cope with biallelic non-SNPs
  # we catch an error as one of the alleles is not a valid allele
  expect_error(vcf_vcfr_gt <- gen_tibble(
    x = vcf_path,
    parser = "vcfR",
    backingfile = tempfile(),
    quiet = TRUE
  ), "valid alleles are ")
})
