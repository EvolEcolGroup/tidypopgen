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

test_that("augment and augment_loci work with a single component", {
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
  augmented_gt <- test_pca %>% augment(data = test_gt, k = 1)
  expect_true(ncol(augmented_gt) == ncol(test_gt) + 1 + 1)
  expect_equal(augmented_gt$.fittedPC1, test_pca$u[, 1] * test_pca$d[1],
    ignore_attr = TRUE
  )
  augmented_gt <- test_pca %>% augment(k = 1)
  expect_true(ncol(augmented_gt) == 1 + 1)
  augmented_loci <- test_pca %>% augment_loci(data = test_gt, k = 1)
  expect_true(ncol(augmented_loci) == ncol(show_loci(test_gt)) + 1)
  augmented_loci <- test_pca %>% augment_loci(k = 1)
  expect_true(ncol(augmented_loci) == 1 + 1)
})

test_that("augment_loci works with gt_pca_autoSVD", {
  set.seed(123)
  n_ind <- 100
  n_loci <- 900
  pops <- rep(1:4, each = n_ind / 4)
  freqs <- matrix(stats::runif(4 * n_loci, 0.05, 0.95), nrow = 4)
  geno <- t(sapply(seq_len(n_ind), function(i) {
    stats::rbinom(n_loci, 2, freqs[pops[i], ])
  }))
  loci <- data.frame(
    name = paste0("rs", seq_len(n_loci)),
    chromosome = rep(1:3, each = n_loci / 3),
    position = rep(seq_len(n_loci / 3) * 1000, 3),
    genetic_dist = 0,
    allele_ref = "A",
    allele_alt = "T"
  )
  meta <- data.frame(id = paste0("ind", seq_len(n_ind)), population = pops)
  test_gt <- gen_tibble(
    geno,
    indiv_meta = meta,
    loci = loci,
    backingfile = tempfile(),
    quiet = TRUE
  )
  test_pca <- gt_pca_autoSVD(test_gt, k = 3, roll_size = 7, verbose = FALSE)
  # some loci are removed by the LD pruning
  expect_true(nrow(test_pca$v) < n_loci)
  augmented_loci <- test_pca %>% augment_loci(data = test_gt)
  expect_identical(augmented_loci$name, loci_names(test_gt))
  # loci used in the pca have their loadings, the others are NA
  used <- match(row_names(test_pca$v), augmented_loci$name)
  expect_equal(augmented_loci$.loadingPC1[used], unclass(test_pca$v)[, 1],
    ignore_attr = TRUE
  )
  expect_true(all(is.na(augmented_loci$.loadingPC1[-used])))
  # error if data does not contain the loci used in the pca
  expect_error(
    test_pca %>% augment_loci(data = test_gt %>% select_loci(1:100)),
    "the loci names in 'data' do not correspond"
  )
})

test_that("augment checks the ids of data", {
  bed_file <- system.file("extdata", "example-missing.bed",
    package = "bigsnpr"
  )
  test_gt <- gen_tibble(
    bed_file,
    backingfile = tempfile("test_"),
    quiet = TRUE
  )
  test_gt <- gt_impute_simple(test_gt, method = "mode")
  # use integer ids
  test_gt$id <- seq_len(nrow(test_gt))
  test_pca <- test_gt %>% gt_pca_randomSVD(k = 5)
  augmented_gt <- test_pca %>% augment(data = test_gt)
  expect_identical(augmented_gt$id, test_gt$id)
  # error if the individuals are not in the same order
  expect_error(
    test_pca %>% augment(data = test_gt[rev(seq_len(nrow(test_gt))), ]),
    "the data id column does not correspond to the individuals"
  )
})
