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

test_that("missing genotypes on the first variant do not set ploidy", {
  # B is missing on the first variant, C is missing on all variants
  vcf_path <- write_small_vcf(list(
    c("1", "10", "rs1", "A", "T", "0/1", ".", "."),
    c("1", "20", "rs2", "C", "G", "0/0", "1/1", "."),
    c("1", "30", "rs3", "G", "A", "1/1", "0/1", ".")
  ))
  vcf_path_tetra <- write_small_vcf(list(
    c("1", "10", "rs1", "A", "T", "0/0/1/1", ".", "."),
    c("1", "20", "rs2", "C", "G", "0/0/0/0", "0/1/1/1", "."),
    c("1", "30", "rs3", "G", "A", "1/1/1/1", "0/0/0/1", ".")
  ))
  for (this_parser in c("cpp", "vcfR")) {
    vcf_gt <- gen_tibble(vcf_path,
      parser = this_parser, backingfile = tempfile(), quiet = TRUE
    )
    expect_equal(show_ploidy(vcf_gt), 2)
    expect_equal(indiv_ploidy(vcf_gt), c(2, 2, 2))
    expect_equal(
      show_genotypes(vcf_gt),
      matrix(c(1, NA, NA, 0, 2, NA, 2, 1, NA), nrow = 3)
    )
    vcf_gt_tetra <- gen_tibble(vcf_path_tetra,
      parser = this_parser, backingfile = tempfile(), quiet = TRUE
    )
    expect_equal(show_ploidy(vcf_gt_tetra), 4)
    expect_equal(indiv_ploidy(vcf_gt_tetra), c(4, 4, 4))
    expect_equal(
      show_genotypes(vcf_gt_tetra),
      matrix(c(2, NA, NA, 0, 3, NA, 4, 1, NA), nrow = 3)
    )
  }
})

test_that("invalid vcf loci do not leave backing files behind", {
  vcf_path <- write_small_vcf(list(
    c("1", "10", "rs1", "A", "T", "0/1", "0/0", "1/1"),
    c("1", "10", "rs2", "C", "G", "0/0", "1/1", "0/1"),
    c("1", "30", "rs3", "G", "A", "1/1", "0/1", "0/1")
  ))
  # vcfR has to remove the backing file, which Windows might prevent
  parsers <- if (.Platform$OS.type == "windows") "cpp" else c("cpp", "vcfR")
  for (this_parser in parsers) {
    backingfile <- tempfile()
    expect_error(
      gen_tibble(vcf_path,
        parser = this_parser, backingfile = backingfile, quiet = TRUE
      ),
      "Your data contain duplicated loci"
    )
    expect_false(file.exists(paste0(backingfile, ".bk")))
    # so the backingfile name can be reused without a new version
    expect_warning(
      vcf_gt <- gen_tibble(vcf_path,
        parser = this_parser, backingfile = backingfile, quiet = TRUE,
        allow_duplicates = TRUE
      ),
      "duplicated loci"
    )
    expect_equal(gt_get_file_names(vcf_gt)[2], paste0(backingfile, ".bk"))
  }
})

test_that("vcf with a trailing blank line or CRLF line endings can be read", {
  vcf_path <- write_small_vcf(list(
    c("1", "10", "rs1", "A", "T", "0/1", "0/0", "1/1"),
    c("1", "20", "rs2", "C", "G", "0/0", "1/1", "0/1"),
    c("1", "30", "rs3", "G", "A", "1/1", "0/1", "./.")
  ))
  vcf_lines <- readLines(vcf_path)
  # trailing blank line
  vcf_blank <- tempfile(fileext = ".vcf")
  writeLines(c(vcf_lines, ""), vcf_blank)
  # CRLF line endings
  vcf_crlf <- tempfile(fileext = ".vcf")
  writeBin(charToRaw(paste0(vcf_lines, "\r\n", collapse = "")), vcf_crlf)
  for (this_parser in c("cpp", "vcfR")) {
    vcf_gt <- gen_tibble(vcf_path,
      parser = this_parser, backingfile = tempfile(), quiet = TRUE
    )
    for (this_vcf in c(vcf_blank, vcf_crlf)) {
      this_gt <- gen_tibble(this_vcf,
        parser = this_parser, backingfile = tempfile(), quiet = TRUE
      )
      expect_identical(this_gt$id, c("A", "B", "C"))
      expect_identical(show_genotypes(this_gt), show_genotypes(vcf_gt))
      expect_identical(show_loci(this_gt), show_loci(vcf_gt))
    }
  }
})

test_that("genotypes with more alleles than max_ploidy are caught", {
  # a tetraploid genotype in the middle of a diploid vcf (the last variant is
  # diploid)
  for (tetra_genotype in c("1/1/1/1", "0/0/0/1")) {
    vcf_path <- write_small_vcf(list(
      c("1", "10", "rs1", "A", "T", "0/1", "0/0", "1/1"),
      c("1", "20", "rs2", "C", "G", "0/1", tetra_genotype, "0/1"),
      c("1", "30", "rs3", "G", "A", "1/1", "0/0", "0/1")
    ))
    for (this_parser in c("cpp", "vcfR")) {
      expect_error(
        gen_tibble(vcf_path,
          parser = this_parser, backingfile = tempfile(), quiet = TRUE
        ),
        "a genotype has more than max_ploidy alleles"
      )
    }
  }
})
