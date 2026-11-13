# limit number of threads for tests
data.table::setDTthreads(2)
if (rlang::is_installed("RhpcBLASctl")) {
  RhpcBLASctl::blas_set_num_threads(2)
  RhpcBLASctl::omp_set_num_threads(2)
}

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
  missing_gt <- gt_impute_simple(missing_gt, method = "mode")
  # we have imputed
  expect_true(gt_has_imputed(missing_gt))
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

test_that("backingfile error if individuals do not match the FBM", {
  # remove an individual from missing_gt
  missing_gt_indiv <- missing_gt[-1, ]
  # try to impute
  expect_error(
    gt_impute_simple(missing_gt_indiv, method = "mode"),
    "The individuals in the gen_tibble do not match the rows of the file"
  )
  # duplicate an individual, keeping the same number of rows as the FBM
  missing_gt_dup <- missing_gt[c(1, seq_len(nrow(missing_gt) - 1)), ]
  expect_identical(nrow(missing_gt_dup), nrow(missing_gt))
  expect_error(
    gt_impute_simple(missing_gt_dup, method = "mode"),
    "use gt_update_backingfile\\(\\) to create a new file backing matrix"
  )
  # neither object was imputed
  expect_false(gt_has_imputed(missing_gt_indiv))
  expect_false(gt_has_imputed(missing_gt_dup))
})

test_that("impute the whole FBM if the gen_tibble has a subset of the loci", {
  # each gen_tibble has its own FBM, as imputation modifies the FBM in place
  all_loci_gt <- gen_tibble(
    bed_file,
    backingfile = tempfile("all_loci_"),
    quiet = TRUE
  )
  all_loci_gt <- gt_impute_simple(all_loci_gt, method = "mode")
  gt_set_imputed(all_loci_gt, set = TRUE)
  all_loci_imputed <- show_genotypes(all_loci_gt)

  # keep a subset of the loci, in a different order
  loci_kept <- c(150:101, 1:50)
  subset_gt <- gen_tibble(
    bed_file,
    backingfile = tempfile("subset_loci_"),
    quiet = TRUE
  ) %>%
    select_loci(all_of(loci_kept))
  # there is missing data to impute in the loci that are kept and removed
  expect_true(anyNA(show_genotypes(subset_gt)))
  expect_true(anyNA(show_genotypes(missing_gt)[, -loci_kept]))
  subset_gt <- gt_impute_simple(subset_gt, method = "mode")
  expect_true(gt_has_imputed(subset_gt))
  # the imputed values are not used by default
  expect_false(gt_uses_imputed(subset_gt))
  expect_true(anyNA(show_genotypes(subset_gt)))
  # with the imputed values, the loci are the same as when imputing all loci
  gt_set_imputed(subset_gt, set = TRUE)
  expect_false(anyNA(show_genotypes(subset_gt)))
  expect_identical(
    show_genotypes(subset_gt),
    all_loci_imputed[, loci_kept]
  )
  # the loci that are not in the gen_tibble were imputed as well
  subset_fbm <- attr(subset_gt$genotypes, "fbm")
  expect_equal(nrow(subset_fbm), nrow(all_loci_gt))
  expect_false(anyNA(subset_fbm[]))
  expect_identical(unname(subset_fbm[]), unname(all_loci_imputed))
  # functions that need imputed data now work on the subset
  expect_true(inherits(subset_gt %>% gt_pca_partialSVD(), "gt_pca"))
})

test_that("impute a gen_tibble with all individuals in a different order", {
  all_indiv_gt <- gen_tibble(
    bed_file,
    backingfile = tempfile("all_indiv_"),
    quiet = TRUE
  )
  imputed_gt <- gt_impute_simple(all_indiv_gt, method = "mode")
  gt_set_imputed(imputed_gt, set = TRUE)
  imputed_genotypes <- show_genotypes(imputed_gt)

  reordered_gt <- gen_tibble(
    bed_file,
    backingfile = tempfile("reordered_indiv_"),
    quiet = TRUE
  )
  new_order <- rev(seq_len(nrow(reordered_gt)))
  reordered_gt <- reordered_gt[new_order, ] %>%
    select_loci(1:200)
  reordered_gt <- gt_impute_simple(reordered_gt, method = "mode")
  gt_set_imputed(reordered_gt, set = TRUE)
  expect_identical(
    show_genotypes(reordered_gt),
    imputed_genotypes[new_order, 1:200]
  )
})


test_that("error imputing an already imputed set", {
  # impute
  missing_gt_imputed <- gt_impute_simple(missing_gt, method = "mode")
  expect_equal(attr(missing_gt_imputed$genotypes, "imputed"), "simple")
  # try to impute again
  expect_error(
    gt_impute_simple(missing_gt_imputed, method = "mode"),
    "object x is already imputed"
  )
  # if the imputation information has been removed (i.e. a corrupted object)
  gt_set_imputed(missing_gt_imputed, TRUE)
  attr(missing_gt_imputed$genotypes, "imputed") <- NULL
  expect_equal(attr(missing_gt_imputed$genotypes, "imputed"), NULL)
  expect_false(gt_has_imputed(missing_gt_imputed))
  expect_error(
    gt_impute_simple(missing_gt_imputed, method = "mode"),
    "^object x is already imputed, but attr"
  )
})

test_that("gt_impute imputes properly", {
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c", "d", "e", "f"),
    population = c("pop1", "pop1", "pop2", "pop1", "pop1", "pop2")
  )

  test_genotypes <- matrix(
    c(
      0,
      2,
      1,
      1,
      0, #
      0,
      0,
      2,
      0,
      1, #
      2,
      0,
      0,
      1,
      1, #
      1,
      1,
      2,
      2,
      2, #
      0,
      0,
      2,
      1,
      1, #
      NA,
      NA,
      NA,
      NA,
      NA
    ),
    nrow = 6,
    byrow = TRUE
  )

  test_loci <- data.frame(
    name = paste0("rs", 1:5),
    chromosome = c(1, 1, 1, 2, 2),
    position = c(3, 65, 343, 23, 456),
    genetic_dist = as.integer(rep(0, 5)),
    allele_ref = c("A", "T", "C", "G", "C"),
    allele_alt = c("T", "C", NA, "C", "G")
  )

  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  )

  # test errors on non-imputed set
  expect_error(
    gt_uses_imputed(test_gt),
    "this dataset does not have any imputed values"
  )
  expect_error(
    gt_set_imputed(test_gt, TRUE),
    "this dataset does not have imputed values"
  )

  # impute method = 'mode'
  imputed_gt_mode <- gt_impute_simple(test_gt, method = "mode")

  # check imputation
  expect_false(gt_has_imputed(test_gt))
  expect_true(gt_has_imputed(imputed_gt_mode))

  # set imputation
  gt_set_imputed(imputed_gt_mode, TRUE)
  expect_false(any(is.na(show_genotypes(imputed_gt_mode))))

  # test error trying to impute an already imputed set
  expect_error(
    gt_impute_simple(imputed_gt_mode),
    "object x is already imputed"
  )

  # Check imputed 'mode' method
  mode_function <- function(x) {
    unique_x <- unique(x)
    return(unique_x[which.max(tabulate(match(x, unique_x)))])
  }

  modes <- apply(test_genotypes, 2, mode_function)
  expect_true(all(show_genotypes(imputed_gt_mode)[6, ] == modes))

  # impute method = 'mean0'
  imputed_gt_mean0 <- gt_impute_simple(test_gt, method = "mean0")

  # set imputation
  gt_set_imputed(imputed_gt_mean0, TRUE)
  expect_false(any(is.na(show_genotypes(imputed_gt_mean0))))

  # check imputed 'mean0' method
  means <- round(colMeans(test_genotypes, na.rm = TRUE), digit = 0)
  expect_true(all(means == show_genotypes(imputed_gt_mean0)[6, ]))

  # impute method = 'random'
  imputed_gt_random <- gt_impute_simple(test_gt, method = "random")

  # set imputation
  gt_set_imputed(imputed_gt_random, TRUE)
  expect_false(any(is.na(show_genotypes(imputed_gt_random))))
})

test_that("n_cores can be set", {
  expect_true(getOption("bigstatsr.check.parallel.blas"))
  bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
  missing_gt <- gen_tibble(
    bed_file,
    backingfile = tempfile("missing_"),
    quiet = TRUE
  )
  one_core <- gt_impute_simple(missing_gt, method = "mode", n_cores = 1)
  two_core <- gt_impute_simple(missing_gt, method = "mode", n_cores = 2)
  expect_true(getOption("bigstatsr.check.parallel.blas"))

  # test parallel blas is true on exit if function errors
  expect_error(
    gt_impute_simple("blah", n_cores = 2),
    "operator is invalid for atomic vectors"
  )
  expect_true(getOption("bigstatsr.check.parallel.blas"))

  # the previous value of the option is restored on exit
  withr::local_options(bigstatsr.check.parallel.blas = FALSE)
  two_core <- gt_impute_simple(missing_gt, method = "mode", n_cores = 2)
  expect_false(getOption("bigstatsr.check.parallel.blas"))
})

test_that("imputed gt retains status after gt_update_backingfile", {
  missing_gt_imputed <- gt_impute_simple(missing_gt, method = "mode")
  # check imputed
  expect_true(gt_has_imputed(missing_gt_imputed))
  # subset by loci
  missing_gt_imputed <- missing_gt_imputed %>%
    select_loci(c(1:100))
  before_update <- missing_gt_imputed
  # update backingfile
  missing_gt_imputed <- gt_update_backingfile(
    missing_gt_imputed,
    tempfile("missing_imputed_"),
    quiet = TRUE
  )
  # check if still imputed
  gt_has_imputed(missing_gt_imputed)
  # the raw values are copied, so genotypes are the same with and without
  # using the imputed values
  expect_identical(
    show_genotypes(missing_gt_imputed),
    show_genotypes(before_update)
  )
  gt_set_imputed(missing_gt_imputed, set = TRUE)
  gt_set_imputed(before_update, set = TRUE)
  expect_identical(
    show_genotypes(missing_gt_imputed),
    show_genotypes(before_update)
  )
  gt_set_imputed(missing_gt_imputed, set = FALSE)
  gt_set_imputed(before_update, set = FALSE)
  # test we can use a function that requires imputation
  expect_true(inherits(missing_gt_imputed %>% gt_pca_partialSVD(), "gt_pca"))


  missing_gt_imputed <- gt_impute_simple(missing_gt, method = "mode")
  # check imputed
  gt_has_imputed(missing_gt_imputed)
  # subset by loci
  missing_gt_imputed <- missing_gt_imputed %>%
    filter(id %in% paste0("ind_", 1:100))
  # update backingfile
  missing_gt_imputed <- gt_update_backingfile(
    missing_gt_imputed,
    tempfile("missing_imputed_"),
    quiet = TRUE
  )
  # check if still imputed
  gt_has_imputed(missing_gt_imputed)
  # test we can use a function that requires imputation
  expect_true(inherits(missing_gt_imputed %>% gt_pca_partialSVD(), "gt_pca"))
})
