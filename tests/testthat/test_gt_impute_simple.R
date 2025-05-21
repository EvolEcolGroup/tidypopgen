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

test_that("backingfile error", {
  # remove an individual from missing_gt
  missing_gt <- missing_gt[-1, ]
  # try to impute
  expect_error(
    missing_gt <- gt_impute_simple(missing_gt, method = "mode"),
    "The number of individuals in the gen_tibble does not match "
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


test_that("imputing subsets", {
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
      NA #
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

  # create a subset of the original tibble
  test_sub <- test_gt[3:6, ]
  test_sub <- test_sub %>% select_loci(c(3:5))

  # impute method = 'mode'
  expect_error(
    imputed_test_sub <- gt_impute_simple(test_sub, method = "mode"),
    "The number of individuals in the gen_tibble does not match"
  )

  # update backingfile
  suppressMessages(test_sub <- gt_update_backingfile(test_sub))

  # try again
  imputed_test_sub <- gt_impute_simple(test_sub, method = "mode")

  # full gen_tibble does not have imputed
  expect_false(gt_has_imputed(test_gt))

  # but the subset does
  expect_true(gt_has_imputed(imputed_test_sub))

  # set imputation
  gt_set_imputed(imputed_test_sub, TRUE)

  # no NA's in the subset
  expect_false(any(is.na(show_genotypes(imputed_test_sub))))

  # but retains NA's in the full gen_tibble
  expect_true(any(is.na(show_genotypes(test_gt))))

  # But it is not possible to have two imputation methods on two subsets of
  # the same data

  # the rest of the individuals and SNPs: 1 and 2
  test_remaining <- test_gt[1:2, ]
  test_remaining <- test_sub %>% select_loci(c(1, 2))

  # impute method = 'mean0'
  impute_remaining <- gt_impute_simple(test_remaining, method = "mean0")
  expect_true(gt_has_imputed(impute_remaining))
  gt_set_imputed(impute_remaining, TRUE)

  # imputes correctly
  expect_false(any(is.na(show_genotypes(impute_remaining))))

  # if we return to our first subset: ind 3:6 and loci 3:5
  expect_false(any(is.na(show_genotypes(imputed_test_sub))))

  mode_function <- function(x) {
    unique_x <- unique(x)
    return(unique_x[which.max(tabulate(match(x, unique_x)))])
  }
  test_sub_genotypes <- test_genotypes[3:6, 3:5]

  modes <- apply(test_sub_genotypes, 2, mode_function)
  means <- round(colMeans(test_sub_genotypes, na.rm = TRUE), digit = 0)

  # expect that the method is overwritten
  expect_false(all(show_genotypes(imputed_test_sub)[4, ] == modes))

  # the original subset now contains imputed means rather than imputed modes
  expect_true(all(show_genotypes(imputed_test_sub)[4, ] == means))
  # we have only changed the other subset (impute_remaining),
  # but imputed_test_sub changes too
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
})
