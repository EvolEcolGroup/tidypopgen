# limit number of threads for tests
data.table::setDTthreads(2)
if (rlang::is_installed("RhpcBLASctl")) {
  RhpcBLASctl::blas_set_num_threads(2)
  RhpcBLASctl::omp_set_num_threads(2)
}

skip_if_not_installed("vcfR")

# write a small vcf; `variants` is a list of character vectors, each with
# chromosome, position, id, ref, alt and the genotypes of samples A, B and C
write_small_vcf <- function(variants, vcf_path = tempfile(fileext = ".vcf")) {
  header <- c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    paste(
      "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
      "FORMAT", "A", "B", "C",
      sep = "\t"
    )
  )
  body <- vapply(variants, function(x) {
    paste(c(x[1:5], ".", "PASS", ".", "GT", x[-(1:5)]), collapse = "\t")
  }, character(1))
  writeLines(c(header, body), vcf_path)
  vcf_path
}

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

  # the vcfR parser also drops biallelic non-SNPs (e.g. indels)
  vcf_vcfr_gt <- gen_tibble(
    x = vcf_path,
    parser = "vcfR",
    backingfile = tempfile(),
    quiet = TRUE
  )
  expect_identical(show_loci(vcf_vcfr_gt), show_loci(vcf_cpp_gt))
  expect_identical(show_genotypes(vcf_vcfr_gt), show_genotypes(vcf_cpp_gt))
})

test_that("partially missing vcf genotypes are coded as missing", {
  vcf_path <- write_small_vcf(list(
    c("1", "10", "rs1", "A", "T", "0/1", "0/0", "1/1"),
    c("1", "20", "rs2", "C", "G", "1/.", "0/.", "./."),
    c("1", "30", "rs3", "G", "A", "1|.", "./1", "0|1")
  ))
  expected <- matrix(c(1, 0, 2, NA, NA, NA, NA, NA, 1), nrow = 3)
  for (this_parser in c("cpp", "vcfR")) {
    vcf_gt <- gen_tibble(vcf_path,
      parser = this_parser, backingfile = tempfile(), quiet = TRUE
    )
    expect_equal(show_genotypes(vcf_gt), expected)
  }
})
