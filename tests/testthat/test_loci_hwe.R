test_that("loci_hwe produces the same output as plink --hardy midp", {
  # Create gentibble for our data
  bed_path <- system.file(
    "extdata/related/families.bed",
    package = "tidypopgen"
  )
  families_bigsnp_path <- bigsnpr::snp_readBed(
    bed_path,
    backingfile = tempfile()
  )
  families <- gen_tibble(
    families_bigsnp_path,
    quiet = TRUE,
    valid_alleles = c("1", "2")
  )

  # Read in plink results (first 10 snps)
  plink_hwe <- read.table(
    system.file(
      "extdata/related/families_hwe_midp.hwe",
      package = "tidypopgen"
    ),
    header = TRUE
  )

  # Calculate hwe in tidypopgen
  tidy_hwe <- families %>%
    select_loci(c(1:10)) %>%
    loci_hwe(mid_p = TRUE)

  # Compare
  result <- all.equal(plink_hwe$P, tidy_hwe, tolerance = 0.0001)

  # Check results are the same to 4 decimals
  expect_true(result)
})


test_that("loci_hwe mid_p = FALSE produces the same output as plink --hardy ", {
  # Create gentibble for our data
  bed_path <- system.file(
    "extdata/related/families.bed",
    package = "tidypopgen"
  )
  families_bigsnp_path <- bigsnpr::snp_readBed(
    bed_path,
    backingfile = tempfile()
  )
  families <- gen_tibble(
    families_bigsnp_path,
    quiet = TRUE,
    valid_alleles = c("1", "2")
  )

  # Read in plink results (first 10 snps)
  plink_hwe <- read.table(
    system.file("extdata/related/families_hwe.hwe", package = "tidypopgen"),
    header = TRUE
  )

  # Calculate hwe in tidypopgen
  tidy_hwe <- families %>%
    select_loci(c(1:10)) %>%
    loci_hwe(mid_p = FALSE)

  # Compare
  result <- all.equal(plink_hwe$P, tidy_hwe, tolerance = 0.0001)

  # Check results are the same to 4 decimals
  expect_true(result)
})
