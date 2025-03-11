test_that("augment_loci adds to loci_table", {
  bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
  missing_gt <- gen_tibble(
    bed_file,
    backingfile = tempfile("missing_"),
    quiet = TRUE
  )
  missing_gt <- gt_impute_simple(missing_gt, method = "mode")
  missing_pca <- missing_gt %>% gt_pca_partialSVD()

  # Add loadings to loci table
  missing_pca_load <- augment_loci(missing_pca, data = missing_gt)
  gt_colnames <- colnames(show_loci(missing_gt))
  expect_true(all(
    colnames(missing_pca_load) ==
      c(
        gt_colnames,
        paste0(".loadingPC", c(1:10))
      )
  ))

  # Try assigning loadings to object with fewer loci
  missing_gt_rm <- missing_gt %>% select_loci(c(1:400))
  expect_error(
    augment_loci(missing_pca, data = missing_gt_rm),
    paste(
      "the loci names in 'data' do not correspond to",
      "the loci in the pca object "
    )
  )
})
