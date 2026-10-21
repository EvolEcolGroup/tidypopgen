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

test_that("count_vcf_individuals counts samples with missing genotypes", {
  vcf_path <- tempfile(fileext = ".vcf")
  writeLines(
    c(
      "##fileformat=VCFv4.2",
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
      paste(
        "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
        "FORMAT", "A", "B", "C",
        sep = "\t"
      ),
      paste("1", "10", "rs1", "A", "T", ".", "PASS", ".", "GT",
        "0/1", "0/0", "1/1",
        sep = "\t"
      ),
      paste("1", "20", "rs2", "C", "G", ".", "PASS", ".", "GT",
        "0/0", ".", "0/1",
        sep = "\t"
      ),
      paste("1", "30", "rs3", "G", "A", ".", "PASS", ".", "GT",
        ".", "0/0", "0/1",
        sep = "\t"
      )
    ),
    vcf_path
  )
  expect_equal(count_vcf_individuals(vcf_path), 3)
  # and the vcfR parser can now read the file
  vcfr_gt <- gen_tibble(vcf_path,
    parser = "vcfR", backingfile = tempfile(), quiet = TRUE
  )
  cpp_gt <- gen_tibble(vcf_path,
    parser = "cpp", backingfile = tempfile(), quiet = TRUE
  )
  expect_equal(vcfr_gt$id, c("A", "B", "C"))
  expect_identical(show_genotypes(vcfr_gt), show_genotypes(cpp_gt))
  expect_equal(
    show_genotypes(vcfr_gt),
    matrix(c(1, 0, 2, 0, NA, 1, NA, 0, 1), nrow = 3)
  )
})
