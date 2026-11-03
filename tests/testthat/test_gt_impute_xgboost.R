# limit number of threads for tests
data.table::setDTthreads(2)
if (rlang::is_installed("RhpcBLASctl")) {
  RhpcBLASctl::blas_set_num_threads(2)
  RhpcBLASctl::omp_set_num_threads(2)
}

skip_if_not_installed("xgboost")

bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
missing_gt <- gen_tibble(
  bed_file,
  backingfile = tempfile("missing_"),
  quiet = TRUE
)

test_that("impute and use the imputation", {
  # we get errors because of missing values
  expect_error(
    missing_gt %>% gt_pca_partialSVD(),
    "You can't have missing values"
  )
  expect_false(gt_has_imputed(missing_gt))
  expect_error(
    gt_uses_imputed(missing_gt),
    "this dataset does not have any imputed"
  )
  # now impute
  missing_gt <- gt_impute_xgboost(missing_gt, seed = 1)
  # we have imputed
  expect_true(gt_has_imputed(missing_gt))
  expect_equal(attr(missing_gt$genotypes, "imputed"), "xgboost")
  # but don't use it by default
  expect_false(gt_uses_imputed(missing_gt))
  # now we return a pca successfully
  expect_true(inherits(missing_gt %>% gt_pca_partialSVD(), "gt_pca"))
  # simple error message
  expect_error(
    gt_set_imputed(missing_gt),
    "set should be either"
  )
})

test_that("gt_impute_xgboost works with loci out of order", {
  set.seed(42)
  n_ind <- 100
  n_loci <- 80
  # correlated blocks of loci, on two chromosomes of different size
  base_geno <- matrix(stats::rbinom(n_ind * 16, 2, 0.5), nrow = n_ind)
  geno <- base_geno[, rep(1:16, each = 5)]
  flip <- matrix(stats::runif(n_ind * n_loci) < 0.05, nrow = n_ind)
  geno[flip] <- stats::rbinom(sum(flip), 2, 0.5)
  geno[matrix(stats::runif(n_ind * n_loci) < 0.1, nrow = n_ind)] <- NA
  loci <- data.frame(
    name = paste0("rs", seq_len(n_loci)),
    chromosome = c(rep(1L, 60), rep(2L, 20)),
    position = c(seq_len(60), seq_len(20)) * 1000L,
    genetic_dist = 0,
    allele_ref = "A",
    allele_alt = "T"
  )
  meta <- data.frame(id = paste0("ind", seq_len(n_ind)), population = "a")
  ordered_gt <- gen_tibble(geno,
    indiv_meta = meta, loci = loci,
    backingfile = tempfile(), quiet = TRUE
  )
  rev_gt <- gen_tibble(geno,
    indiv_meta = meta, loci = loci,
    backingfile = tempfile(), quiet = TRUE
  )
  # reverse the order of the loci table
  rev_gt <- rev_gt %>% select_loci(rev(seq_len(n_loci)))
  ordered_gt <- gt_impute_xgboost(ordered_gt, seed = 1)
  rev_gt <- gt_impute_xgboost(rev_gt, seed = 1)
  gt_set_imputed(ordered_gt, set = TRUE)
  gt_set_imputed(rev_gt, set = TRUE)
  # the imputed genotypes should be the same
  expect_identical(
    show_genotypes(rev_gt)[, rev(seq_len(n_loci))],
    show_genotypes(ordered_gt)
  )
})

test_that("backingfile error if gt and bk have uneven number of loci/indivs", {
  # remove an individual from missing_gt
  missing_indiv_gt <- missing_gt[-1, ]
  # try to impute
  expect_error(
    missing_indiv_gt <- gt_impute_xgboost(missing_indiv_gt, seed = 1),
    "The number of individuals in the gen_tibble does not match "
  )

  # remove loci from missing_gt
  missing_loci_gt <- missing_gt %>% select_loci(1:5)
  expect_error(
    missing_gt <- gt_impute_xgboost(missing_loci_gt, seed = 1),
    "The number of loci in the gen_tibble does not match "
  )
})

test_that("error imputing an already imputed set", {
  # impute
  missing_gt_imputed <- gt_impute_xgboost(missing_gt, seed = 1)
  expect_equal(attr(missing_gt_imputed$genotypes, "imputed"), "xgboost")
  # try to impute again
  expect_error(
    gt_impute_xgboost(missing_gt_imputed, seed = 1),
    "object x is already imputed"
  )
  # corrupted file
  gt_set_imputed(missing_gt_imputed, TRUE)
  attr(missing_gt_imputed$genotypes, "imputed") <- NULL
  expect_equal(
    attr(missing_gt_imputed$genotypes, "imputed", exact = TRUE),
    NULL
  )
  expect_false(gt_has_imputed(missing_gt_imputed))
  expect_error(
    gt_impute_xgboost(missing_gt_imputed, seed = 1),
    "^object x is already imputed, but attr"
  )
})

test_that("n_cores can be set", {
  bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
  missing_gt <- gen_tibble(
    bed_file,
    backingfile = tempfile("missing_"),
    quiet = TRUE
  )
  expect_true(getOption("bigstatsr.check.parallel.blas"))
  one_core <- gt_impute_xgboost(missing_gt, n_cores = 1, seed = 1)
  two_core <- gt_impute_xgboost(missing_gt, n_cores = 2, seed = 1)
  expect_true(getOption("bigstatsr.check.parallel.blas"))

  # test parallel blas is true on exit if function errors
  expect_error(
    gt_impute_xgboost("blah", n_cores = 2),
    "operator is invalid for atomic vectors"
  )
  expect_true(getOption("bigstatsr.check.parallel.blas"))
})

test_that("append_error correct dimensions", {
  bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
  missing_gt <- gen_tibble(
    bed_file,
    backingfile = tempfile("missing_"),
    quiet = TRUE
  )
  missing_gt <- gt_impute_xgboost(missing_gt, seed = 1, append_error = TRUE)
  expect_true(gt_has_imputed(missing_gt))
  expect_true("imputed_errors" %in% names(attributes(missing_gt$genotypes)))
  expect_equal(dim(attr(missing_gt$genotypes, "imputed_errors")), c(2, 500))
})
