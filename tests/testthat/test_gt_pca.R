testthat::test_that("gt_pca is equivalent to adegenet glPca",{
library(adegenet)
  # simulate object
  set.seed(123)
  x <- glSim(50,4e3, 50, ploidy=2)
  # fit pca in adegenet
  pca_gl <- glPca(x, nf=2)
  # convert object to gen_tibble
  x_gt <- as_gen_tibble(x)
  # use gt_pca
  pca_gt <- gt_pca(x_gt, nf=2)
  expect_true(all.equal(pca_gl$scores,pca_gt$scores, check.attributes = FALSE))
  # now use C to compute the dot product
  pca_gt_c <- gt_pca(x_gt, nf=2, useC = TRUE)
  expect_identical(pca_gt$scores, pca_gt_c$scores)
})
