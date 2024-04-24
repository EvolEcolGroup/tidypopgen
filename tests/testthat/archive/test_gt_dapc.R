testthat::test_that("gt_pca is equivalent to adegenet glPca",{
  library(adegenet)
  # simulate object
  set.seed(123)
  x <- glSim(50,4e3, 50, ploidy=2)
  x_gt <- as_gen_tibble(x)
  # fit gt_pca
  test_pca <- gt_pca(x_gt, center = TRUE)
  test_clusters <- gt_cluster_pca(test_pca)
  # in the assignments for k=3, the maximum group id should be 3
  expect_true(max(test_clusters$clusters$groups[[3]])==3)
  test_clusters <- gt_cluster_pca_best_k(test_clusters, quiet=TRUE)
  # check that we now have the info
  expect_true(!is.null(test_clusters$best_k))
  # and now dapc
  test_dapc <- gt_dapc(test_clusters)
  # simple test
  expect_true(nlevels(test_dapc$grp)==test_clusters$best_k)
  # we should test that we get the expected results when we set pop manually
})
