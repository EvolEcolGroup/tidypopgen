options(mc_doScale_quiet = TRUE)

test_that("gt_pca_tidiers", {
  bed_file <- system.file("extdata", "example-missing.bed",
    package = "bigsnpr"
  )
  test_gt <- gen_tibble(
    bed_file,
    backingfile = tempfile("test_"),
    quiet = TRUE
  )

  test_gt <- gt_impute_simple(test_gt, method = "mode")
  test_pca <- test_gt %>% gt_pca_randomSVD(k = 5)
  # get error if more than one matrix selected
  expect_error(
    test_pca %>% tidy(matrix = c("X", "X")),
    "Must select a single matrix to tidy"
  )
  # get eigenvalues
  test_eigenvalues <- test_pca %>% tidy(matrix = "eigenvalues")
  # expect 5 eigenvalues
  expect_true(nrow(test_eigenvalues) == 5)
  # get loadings
  test_loadings <- test_pca %>% tidy(matrix = "loadings")
  # expect loadings to be as many as loci times number of PCs
  expect_true(nrow(show_loci(test_gt)) * 5 == nrow(test_loadings))
  # get scores
  test_scores <- test_pca %>% tidy(matrix = "scores")
  # expect scores to be as many as samples times number of PCs
  expect_true(nrow(test_gt) * 5 == nrow(test_scores))
  #######################
  # now augment
  augmented_gt <- test_pca %>% augment(data = test_gt)
  # expect augmented to have same number of rows as original gt
  expect_true(nrow(augmented_gt) == nrow(test_gt))
  # expect augmented to have 5 additional columns plus one for .rownames
  expect_true(ncol(augmented_gt) == ncol(test_gt) + 5 + 1)
  # now augment wihtout data
  augmented_gt <- test_pca %>% augment()
  # expect augmented to have 5 additional columns plus one for .rownames
  expect_true(ncol(augmented_gt) == 5 + 1)
  # now augment with a different number of PCs
  augmented_gt <- test_pca %>% augment(data = test_gt, k = 3)
  # expect augmented to have 3 additional columns plus one for .rownames
  expect_true(ncol(augmented_gt) == ncol(test_gt) + 3 + 1)
  # now augment with a different number of PCs and no data
  augmented_gt <- test_pca %>% augment(k = 3)
  # expect augmented to have 3 additional columns plus one for .rownames
  expect_true(ncol(augmented_gt) == 3 + 1)
  #######################
  # now augment loci
  augmented_loci <- test_pca %>% augment_loci(data = test_gt)
  # expect augmented to have same number of rows as the loci in the original gt
  expect_true(nrow(augmented_loci) == nrow(show_loci(test_gt)))
  # expect augmented to have 5 additional columns
  # compared to show_loci
  expect_true(ncol(augmented_loci) == ncol(show_loci(test_gt)) + 5)
  # now augment loci without data
  augmented_loci <- test_pca %>% augment_loci()
  # expect augmented to have 5 columns plus .rownames
  expect_true(ncol(augmented_loci) == 5 + 1)
  # now augment loci with a different number of PCs
  augmented_loci <- test_pca %>% augment_loci(data = test_gt, k = 3)
  # expect augmented to have 3 additional columns
  # compared to show_loci
  expect_true(ncol(augmented_loci) == ncol(show_loci(test_gt)) + 3)
  # now augment loci with a different number of PCs and no data
  augmented_loci <- test_pca %>% augment_loci(k = 3)
  # expect augmented to have 3 columns plus .rownames
  expect_true(ncol(augmented_loci) == 3 + 1)
})
