testthat::test_that("read_loc_raw correctly reads loc information", {
  # a simple file
  raw_path <- system.file("extdata/pop_b.raw", package = "tidypopgen")
  map_path <- system.file("extdata/pop_b.map", package = "tidypopgen")
  pop_b_gen <- read_plink_raw(file = raw_path, map_file = map_path, quiet = TRUE)
  pop_b_loc_info <- read_plink_raw_loci(raw_path)
  # check that loci names are the same irrespective of how we get them

  expect_identical(show_loci(pop_b_gen) %>% select(name,allele_ref,allele_alt),
                   pop_b_loc_info)
})
