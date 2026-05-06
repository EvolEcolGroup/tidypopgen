test_that("using integers as names works", {
  vcf_path <-
    system.file("/extdata/anolis/punctatus_t70_s10_n46_filtered.recode.vcf.gz",
      package = "tidypopgen"
    )
  anole_gt <-
    gen_tibble(vcf_path,
      quiet = TRUE, backingfile = tempfile("anolis_"),
      names_as_int = TRUE
    )
  pops_path <- system.file("/extdata/anolis/punctatus_n46_meta.csv",
    package = "tidypopgen"
  )
  pops <- read.csv(pops_path)
  anole_gt <- anole_gt %>% left_join(pops, by = "id")
  # check that names are integers
  expect_true(is.integer(show_loci(anole_gt)$name))
  anole_gt <- gt_add_sf(anole_gt, c("longitude", "latitude"))
  expect_true(is.integer(show_loci(anole_gt)$name))
  # test a pca
  anole_gt <- gt_impute_simple(anole_gt, method = "mode")
  anole_pca <- anole_gt %>% gt_pca_randomSVD(k = 30)
  expect_true(is.integer(row_names(anole_pca$v)))
  anole_pca <- anole_gt %>% gt_pca_partialSVD(k = 30)
  expect_true(is.integer(row_names(anole_pca$v)))
  # tidy the pca
  loadings_pca <- tidy(anole_pca, matrix = "loadings")
  # check that we have a an integer column for "column"
  expect_true(is.integer(loadings_pca$column))
  # TODO check what happens if we have snps not in perfect order 1 to X (i.e.
  # some missing numbers) we might need to implement a method for as.data.frame
  # that transfers row_names and the same for as_tibble()
})
