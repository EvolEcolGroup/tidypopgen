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
  gt_set_imputed(test_gt, set = TRUE)
  test_pca <- test_gt %>% gt_pca(k = 15, n_cores = 2)
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
  #######################
  # loadings of several axes are plotted with geom_hex(), which needs hexbin
  expect_s3_class(autoplot(test_dapc, type = "loadings", ld = 1), "ggplot")
  if (rlang::is_installed("hexbin")) {
    expect_s3_class(
      autoplot(test_dapc, type = "loadings", ld = c(1, 2)),
      "patchwork"
    )
  } else {
    rlang::local_interactive(FALSE)
    expect_error(
      autoplot(test_dapc, type = "loadings", ld = c(1, 2)),
      class = "rlib_error_package_not_found"
    )
  }
})

test_that("scores plot error if only 2 populations", {
  bed_file <- system.file("extdata", "example-missing.bed",
    package = "bigsnpr"
  )
  test_gt <- gen_tibble(
    bed_file,
    backingfile = tempfile("test_"),
    quiet = TRUE
  )
  test_gt$population <- c(rep("pop1", 100), rep("pop2", 100))
  test_gt <- test_gt %>% group_by(population)
  test_gt <- gt_impute_simple(test_gt, method = "mode")
  gt_set_imputed(test_gt, set = TRUE)
  test_pca <- test_gt %>% gt_pca(k = 2, n_cores = 2)
  test_dapc <- test_pca %>% gt_dapc(pop = as.factor(test_gt$population))
  expect_error(
    autoplot(test_dapc, type = "scores"),
    "object has only one discriminant axis, cannot plot scores"
  )
})

test_that("gt_cluster_pca works with a single k and a single component", {
  bed_file <- system.file("extdata", "example-missing.bed",
    package = "bigsnpr"
  )
  test_gt <- gen_tibble(
    bed_file,
    backingfile = tempfile("test_"),
    quiet = TRUE
  )
  test_gt <- gt_impute_simple(test_gt, method = "mode")
  gt_set_imputed(test_gt, set = TRUE)
  test_pca <- test_gt %>% gt_pca(k = 5, n_cores = 2)
  # only one k larger than 1
  test_cluster <- gt_cluster_pca(test_pca,
    n_pca = 5, k_clusters = c(1, 2),
    quiet = TRUE
  )
  expect_identical(names(test_cluster$clusters$groups), c("1", "2"))
  expect_identical(
    names(test_cluster$clusters$groups[["1"]]),
    rownames(test_pca$u)
  )
  # only k = 1
  test_cluster <- gt_cluster_pca(test_pca,
    n_pca = 5, k_clusters = 1,
    quiet = TRUE
  )
  expect_identical(names(test_cluster$clusters$groups), "1")
  expect_true(all(test_cluster$clusters$groups[["1"]] == 1))
  # a single component
  test_cluster <- gt_cluster_pca(test_pca,
    n_pca = 1, k_clusters = c(1, 3),
    quiet = TRUE
  )
  expect_length(test_cluster$clusters$WSS, 3)
  expect_identical(
    names(test_cluster$clusters$groups[["3"]]),
    rownames(test_pca$u)
  )
})

test_that("best_k is a value of k when k_clusters does not start at 1", {
  bed_file <- system.file("extdata", "example-missing.bed",
    package = "bigsnpr"
  )
  test_gt <- gen_tibble(
    bed_file,
    backingfile = tempfile("test_"),
    quiet = TRUE
  )
  test_gt <- gt_impute_simple(test_gt, method = "mode")
  gt_set_imputed(test_gt, set = TRUE)
  test_pca <- test_gt %>% gt_pca(k = 10, n_cores = 2)
  test_cluster <- gt_cluster_pca(test_pca,
    n_pca = 10, k_clusters = c(3, 8),
    quiet = TRUE
  )
  for (this_criterion in c("min", "diffNgroup")) {
    test_cluster_best <- gt_cluster_pca_best_k(test_cluster,
      stat = "WSS", criterion = this_criterion, quiet = TRUE
    )
    expect_true(test_cluster_best$best_k %in% test_cluster$clusters$k)
  }
  test_cluster_best <- gt_cluster_pca_best_k(test_cluster,
    stat = "WSS", criterion = "min", quiet = TRUE
  )
  # the WSS decreases with k, so the minimum is at the largest k
  expect_equal(
    test_cluster_best$best_k,
    test_cluster$clusters$k[which.min(test_cluster$clusters$WSS)]
  )
  expect_message(
    gt_cluster_pca_best_k(test_cluster, stat = "WSS", criterion = "min"),
    paste0(test_cluster_best$best_k, " clusters")
  )
  # gt_dapc uses the groups for best_k
  test_dapc <- gt_dapc(test_cluster_best)
  expect_equal(nlevels(test_dapc$grp), test_cluster_best$best_k)
  # best_k can be overridden by the user
  test_cluster_best$best_k <- 4
  test_dapc <- gt_dapc(test_cluster_best)
  expect_equal(nlevels(test_dapc$grp), 4)
  # pop can be given as the value of k
  test_dapc <- gt_dapc(test_cluster, pop = 5)
  expect_equal(nlevels(test_dapc$grp), 5)
})

test_that("gt_dapc var is the proportion of variance of the retained PCs", {
  bed_file <- system.file("extdata", "example-missing.bed",
    package = "bigsnpr"
  )
  test_gt <- gen_tibble(
    bed_file,
    backingfile = tempfile("test_"),
    quiet = TRUE
  )
  test_gt <- gt_impute_simple(test_gt, method = "mode")
  gt_set_imputed(test_gt, set = TRUE)
  test_pca <- test_gt %>% gt_pca(k = 10, n_cores = 2)
  pops <- as.factor(rep(c("pop1", "pop2", "pop3", "pop4"), each = 50))
  test_dapc <- gt_dapc(test_pca, pop = pops, n_pca = 6)
  # the eigenvalues are the squared singular values
  expect_equal(test_dapc$var, sum(test_pca$d[1:6]^2) / sum(test_pca$d^2))
})

test_that("gt_dapc uses k-1 PCs by default", {
  bed_file <- system.file("extdata", "example-missing.bed",
    package = "bigsnpr"
  )
  test_gt <- gen_tibble(
    bed_file,
    backingfile = tempfile("test_"),
    quiet = TRUE
  )
  test_gt <- gt_impute_simple(test_gt, method = "mode")
  gt_set_imputed(test_gt, set = TRUE)
  test_pca <- test_gt %>% gt_pca(k = 5, n_cores = 2)
  # as many populations as PCs
  pops <- as.factor(rep(paste0("pop", 1:5), each = 40))
  test_dapc <- gt_dapc(test_pca, pop = pops)
  expect_equal(test_dapc$n.pca, 4)
  # fewer populations than PCs
  pops <- as.factor(rep(paste0("pop", 1:3), length.out = 200))
  test_dapc <- gt_dapc(test_pca, pop = pops)
  expect_equal(test_dapc$n.pca, 2)
  # more populations than PCs
  pops <- as.factor(rep(paste0("pop", 1:8), each = 25))
  test_dapc <- gt_dapc(test_pca, pop = pops)
  expect_equal(test_dapc$n.pca, 5)
})

test_that("gt_dapc gives an informative error for an invalid pop", {
  bed_file <- system.file("extdata", "example-missing.bed",
    package = "bigsnpr"
  )
  test_gt <- gen_tibble(
    bed_file,
    backingfile = tempfile("test_"),
    quiet = TRUE
  )
  test_gt <- gt_impute_simple(test_gt, method = "mode")
  gt_set_imputed(test_gt, set = TRUE)
  test_pca <- test_gt %>% gt_pca(k = 5, n_cores = 2)
  pops <- c(rep("pop1", 100), rep("pop2", 100))
  # a character vector rather than a factor
  expect_error(
    gt_dapc(test_pca, pop = pops),
    "'pop' should be a factor, or a value of k"
  )
  # a value of k, but x is not a gt_cluster_pca
  expect_error(
    gt_dapc(test_pca, pop = 3),
    "'pop' should be a factor, or a value of k"
  )
  # a value of k that was not explored
  test_cluster <- gt_cluster_pca(test_pca,
    n_pca = 5, k_clusters = c(1, 4),
    quiet = TRUE
  )
  expect_error(
    gt_dapc(test_cluster, pop = 10),
    "'pop' should be a factor, or a value of k"
  )
})

test_that("augment.gt_dapc works with one discriminant axis and checks ids", {
  bed_file <- system.file("extdata", "example-missing.bed",
    package = "bigsnpr"
  )
  test_gt <- gen_tibble(
    bed_file,
    backingfile = tempfile("test_"),
    quiet = TRUE
  )
  test_gt <- gt_impute_simple(test_gt, method = "mode")
  gt_set_imputed(test_gt, set = TRUE)
  # use integer ids
  test_gt$id <- seq_len(nrow(test_gt))
  test_pca <- test_gt %>% gt_pca(k = 5, n_cores = 2)
  pops <- as.factor(c(rep("pop1", 100), rep("pop2", 100)))
  test_dapc <- test_pca %>% gt_dapc(pop = pops, n_pca = 3)
  expect_equal(test_dapc$n.da, 1)
  augmented_gt <- test_dapc %>% augment(data = test_gt)
  expect_true(ncol(augmented_gt) == ncol(test_gt) + 1 + 1)
  expect_equal(augmented_gt$.fittedLD1, test_dapc$ind.coord[, 1],
    ignore_attr = TRUE
  )
  augmented_gt <- test_dapc %>% augment()
  expect_true(ncol(augmented_gt) == 1 + 1)
  # error if the individuals are not in the same order
  expect_error(
    test_dapc %>% augment(data = test_gt[rev(seq_len(nrow(test_gt))), ]),
    "the data id column does not correspond to the individuals"
  )
})

test_that("gt_dapc with a PCA of missing genotypes", {
  bed_file <- system.file("extdata", "example-missing.bed",
    package = "bigsnpr"
  )
  test_gt <- gen_tibble(
    bed_file,
    backingfile = tempfile("test_"),
    quiet = TRUE
  )
  pops <- as.factor(c(rep("pop1", 100), rep("pop2", 100)))
  # a PCA that compared individuals over the loci genotyped in both
  pairwise_pca <- test_gt %>% gt_pca(k = 5, n_cores = 2)
  # the loadings of the loci need the genotypes
  expect_error(
    gt_dapc(pairwise_pca, pop = pops, n_pca = 3),
    "the loadings of the loci need the genotypes"
  )
  # without the loadings of the loci, the genotypes are not needed
  no_loadings <- gt_dapc(pairwise_pca,
    pop = pops, n_pca = 3, loadings_by_locus = FALSE
  )
  expect_null(no_loadings$var.load)
  # with them, the loadings of the loci are their regression coefficients on
  # the scores, multiplied by the loadings of the discriminant functions
  pairwise_dapc <- gt_dapc(pairwise_pca,
    pop = pops, n_pca = 3, data = test_gt, n_cores = 2
  )
  expect_identical(pairwise_dapc$ind.coord, no_loadings$ind.coord)
  lsq <- loci_lsq_loadings(pairwise_pca, test_gt, n_pca = 3, n_cores = 2)
  expect_equal(
    as.vector(unclass(pairwise_dapc$var.load)),
    as.vector(lsq %*% pairwise_dapc$loadings)
  )
  expect_identical(row_names(pairwise_dapc$var.load), loci_names(test_gt))
  # also through the clusters made from the PCA
  pairwise_clusters <- gt_cluster_pca(pairwise_pca,
    n_pca = 5, k_clusters = c(1, 3), quiet = TRUE
  )
  expect_s3_class(
    gt_dapc(pairwise_clusters, pop = 2, data = test_gt),
    "gt_dapc"
  )
  # data must be the gen_tibble of the PCA
  expect_error(
    gt_dapc(pairwise_pca, pop = pops, data = test_gt[200:1, ]),
    "the individuals in 'data' do not correspond to those in the PCA"
  )
  expect_error(
    gt_dapc(pairwise_pca, pop = pops, data = select_loci(test_gt, 1:10)),
    "loci used in the PCA are not present in 'data'"
  )
  # with a PCA of imputed genotypes, the loadings of the loci are those of
  # the PCA, multiplied by the loadings of the discriminant functions (as
  # before), and do not need the genotypes
  test_gt <- gt_impute_simple(test_gt, method = "mode")
  gt_set_imputed(test_gt, set = TRUE)
  imputed_pca <- test_gt %>% gt_pca(k = 5, n_cores = 2)
  imputed_dapc <- gt_dapc(imputed_pca, pop = pops, n_pca = 3)
  expect_equal(
    as.vector(unclass(imputed_dapc$var.load)),
    as.vector(unclass(imputed_pca$v)[, 1:3] %*% imputed_dapc$loadings)
  )
})

test_that("loadings of loci by least squares with missing genotypes", {
  set.seed(1)
  sim_gt <- gt_simulate_bn(
    n_indiv = c(20, 20, 20), n_loci = 400, fst = 0.1,
    missing = rep(seq(0, 0.6, length.out = 20), 3)
  )
  sim_gt <- select_loci_if(sim_gt, loci_maf(genotypes) > 0)
  sim_pca <- gt_pca(sim_gt, k = 4, n_cores = 1)
  n_pca <- 3
  lsq <- loci_lsq_loadings(sim_pca, sim_gt, n_pca = n_pca, n_cores = 2)
  # the reference: the regression of the standardised genotypes of each
  # locus on an intercept and the scores, over the individuals genotyped at
  # that locus
  geno <- show_genotypes(sim_gt)
  z <- sweep(sweep(geno, 2, sim_pca$center), 2, sim_pca$scale, "/")
  scores <- sweep(sim_pca$u[, 1:n_pca], 2, sim_pca$d[1:n_pca], "*")
  ref <- t(sapply(seq_len(ncol(z)), function(l) {
    obs <- !is.na(z[, l])
    stats::lm.fit(cbind(1, scores[obs, ]), z[obs, l])$coefficients[-1]
  }))
  expect_equal(unname(lsq), unname(ref), tolerance = 1e-8)
  # the BLAS and the internal engine give the same coefficients (with several
  # blocks of loci)
  args <- list(
    BM = .gt_get_fbm(sim_gt), rowInd = .gt_fbm_rows(sim_gt),
    colInd = .gt_fbm_cols(sim_gt), center = sim_pca$center,
    scale = sim_pca$scale, S = unclass(sim_pca$u)[, 1:n_pca],
    block_size = 64L, ncores = 2L
  )
  expect_equal(
    do.call(grm_loci_lsq_cpp, c(args, use_blas = TRUE)),
    do.call(grm_loci_lsq_cpp, c(args, use_blas = FALSE)),
    tolerance = 1e-10
  )

  # with complete genotypes, they are the loadings of the PCA
  complete_gt <- gt_simulate_bn(
    n_indiv = c(20, 20, 20), n_loci = 400, fst = 0.1
  )
  complete_gt <- select_loci_if(complete_gt, loci_maf(genotypes) > 0)
  complete_pca <- gt_pca(complete_gt, k = 4, n_cores = 1)
  expect_equal(
    unname(loci_lsq_loadings(complete_pca, complete_gt,
      n_pca = n_pca, n_cores = 2
    )),
    unname(unclass(complete_pca$v)[, 1:n_pca]),
    tolerance = 1e-6
  )

  # a locus genotyped in fewer individuals than there are coefficients
  # (n_pca + 1) has NA loadings, with a warning
  few_gt <- gen_tibble(
    x = rbind(
      c(2, 1, 0, 1, 1, NA), c(0, 1, 2, 0, 1, NA), c(1, 2, 0, 1, 0, NA),
      c(2, 0, 1, 2, 1, 1), c(0, 2, 1, 0, 2, 0), c(1, 1, 2, 2, 0, NA),
      c(2, 0, 0, 1, 2, NA), c(0, 1, 1, 0, 1, NA)
    ),
    indiv_meta = data.frame(
      id = paste0("ind", 1:8), population = rep(c("a", "b"), each = 4)
    ),
    loci = data.frame(
      name = paste0("rs", 1:6), chromosome = 1, position = 1:6 * 100,
      genetic_dist = 0, allele_ref = "A", allele_alt = "T"
    ),
    backingfile = tempfile(), quiet = TRUE
  )
  few_pca <- gt_pca(few_gt, k = 2, n_cores = 1)
  expect_warning(
    few <- loci_lsq_loadings(few_pca, few_gt, n_pca = 2, n_cores = 1),
    "1 loci are genotyped in too few individuals"
  )
  expect_true(all(is.na(few[6, ])))
  expect_false(anyNA(few[1:5, ]))
})
