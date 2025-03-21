options(mc_doScale_quiet = TRUE)

test_that("gt_dapc_and_tidiers", {
  bed_file <- system.file("extdata", "example-missing.bed",
    package = "bigsnpr"
  )
  test_gt <- gen_tibble(
    bed_file,
    backingfile = tempfile("test_"),
    quiet = TRUE
  )

  test_gt <- gt_impute_simple(test_gt, method = "mode")
  test_pca <- test_gt %>% gt_pca_randomSVD(k = 15)
  ############################
  # unit tests for gt_cluster_pca
  # test error if not a gt object
  expect_error(
    gt_cluster_pca(test_gt),
    "'x' should be a 'gt_pca' object"
  )
  # use all pca
  expect_message(
    test_cluster <- gt_cluster_pca(test_pca),
    "'n.pca' was not set: all 15 components were used"
  )
  # we have a clusters element in the object
  expect_true("clusters" %in% names(test_cluster))
  # TODO we should test that all the sub-element are what we expect

  ############################
  # unit tests for gt_cluster_pca_best_k

  # test error if not a gt_cluster_pca object
  expect_error(
    test_pca %>% gt_cluster_pca_best_k(stat = "X", criterion = "min"),
    "'x' should be a 'gt_cluster_pca' object"
  )
  expect_message(
    test_cluster_best <-
      test_cluster %>%
      gt_cluster_pca_best_k(
        stat = "AIC",
        criterion = "diffNgroup"
      ),
    "Using AIC with criterion"
  )
  # we have a best_k element
  expect_true("best_k" %in% names(test_cluster_best))

  ############################
  # unit tests for gt_dapc
  test_dapc <- test_cluster_best %>% gt_dapc()
  # Test that input validation works
  expect_error(
    gt_dapc("not_a_cluster_object"),
    "'x' should be a 'gt_pca' object"
  )
  # Test that output has expected structure
  expect_true("eig" %in% names(test_dapc))
  expect_true("loadings" %in% names(test_dapc))
  expect_true("ind.coord" %in% names(test_dapc))
  expect_equal(ncol(test_dapc$ind.coord), 5)

  # Test with different n.da parameter
  test_dapc_custom <- test_cluster_best %>% gt_dapc(n_da = 3)
  expect_equal(ncol(test_dapc_custom$ind.coord), 3)


  ############################
  # test the tidiers

  # get error if more than one matrix selected
  expect_error(
    test_dapc %>% tidy(matrix = c("X", "X")),
    "Must select a single matrix to tidy"
  )
  # get eigenvalues
  test_eigenvalues <- test_dapc %>% tidy(matrix = "eigenvalues")
  # expect 5 eigenvalues
  expect_true(nrow(test_eigenvalues) == 5)
  # get loadings
  test_loadings <- test_dapc %>% tidy(matrix = "loadings")
  # expect loadings to be as many as loci times number of PCs
  expect_true(nrow(show_loci(test_gt)) * 5 == nrow(test_loadings))
  # get scores
  test_scores <- test_dapc %>% tidy(matrix = "scores")
  # expect scores to be as many as samples times number of PCs
  expect_true(nrow(test_gt) * 5 == nrow(test_scores))
  #######################
  # now augment
  augmented_gt <- test_dapc %>% augment(data = test_gt)
  # expect augmented to have same number of rows as original gt
  expect_true(nrow(augmented_gt) == nrow(test_gt))
  # expect augmented to have 5 additional columns plus one for .rownames
  expect_true(ncol(augmented_gt) == ncol(test_gt) + 5 + 1)
  # now augment wihtout data
  augmented_gt <- test_dapc %>% augment()
  # expect augmented to have 5 additional columns plus one for .rownames
  expect_true(ncol(augmented_gt) == 5 + 1)
  # now augment with a different number of PCs
  augmented_gt <- test_dapc %>% augment(data = test_gt, k = 3)
  # expect augmented to have 3 additional columns plus one for .rownames
  expect_true(ncol(augmented_gt) == ncol(test_gt) + 3 + 1)
  # now augment with a different number of PCs and no data
  augmented_gt <- test_dapc %>% augment(k = 3)
  # expect augmented to have 3 additional columns plus one for .rownames
  expect_true(ncol(augmented_gt) == 3 + 1)
})
