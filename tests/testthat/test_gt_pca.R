options(mc_doScale_quiet = TRUE)

test_that("fit_gt_pca_and_predict", {
  bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
  missing_gt <- gen_tibble(
    bed_file,
    backingfile = tempfile("missing_"),
    quiet = TRUE
  )
  expect_error(
    missing_gt %>% gt_pca_partialSVD(),
    "You can't have missing"
  )
  missing_gt <- gt_impute_simple(missing_gt, method = "mode")
  missing_pca <- missing_gt %>% gt_pca_partialSVD()
  # check that predicting on the object is the same as predicting from the
  # full dataset without imputation to the center (the data are already imputed)
  expect_true(all.equal(
    predict(missing_pca),
    predict(missing_pca, new_data = missing_gt),
    check.attributes = FALSE
  ))
  # now mismatch the loci table
  missing_gt_edited <- missing_gt
  show_loci(missing_gt_edited)$name[3] <- "blah"
  expect_error(
    predict(missing_pca, new_data = missing_gt_edited),
    "loci used in object"
  )
  missing_gt_edited <- missing_gt
  show_loci(missing_gt_edited)$allele_ref[3] <- "blah"
  expect_error(
    predict(missing_pca, new_data = missing_gt_edited),
    "ref and alt alleles differ"
  )
  # predict when new dataset has extra positions
  missing_gt_sub <- missing_gt %>% select_loci(100:450)
  missing_sub_pca <- missing_gt_sub %>% gt_pca_partialSVD()
  expect_true(all.equal(
    predict(missing_sub_pca),
    predict(missing_sub_pca, new_data = missing_gt, project_method = "none"),
    check.attributes = FALSE
  ))
})

test_that("adjusting roll_size fixes gt_pca_autoSVD problem ", {
  # Generate example_indiv_meta data frame
  individual_ids <- paste0("indiv", 1:200)
  populations <- sample(
    c("pop1", "pop2", "pop3", "pop4"),
    size = 200,
    replace = TRUE
  )
  example_indiv_meta <- data.frame(
    id = individual_ids,
    population = populations
  )

  # Generate example_genotypes matrix (no missingness)
  values <- c(0, 1, 2)
  genotype_matrix <- sample(values, size = 200 * 500, replace = TRUE)
  example_genotypes <- matrix(genotype_matrix, nrow = 200, ncol = 500)

  # Generate example_loci data frame
  loci_names <- paste0("locus", 1:500)
  chromosomes <- c(rep(1, 255), rep(2, 245))
  positions <- seq(from = 1000, by = 1000, length.out = 500)
  allele_refs <- sample(c("A", "T", "C", "G"), size = 500, replace = TRUE)
  allele_alts <- sample(c("A", "T", "C", "G", NA), size = 500, replace = TRUE)

  example_loci <- data.frame(
    name = loci_names,
    chromosome = chromosomes,
    position = positions,
    genetic_dist = rep(0, 500),
    allele_ref = allele_refs,
    allele_alt = allele_alts
  )

  # create gen_tibble
  test <- gen_tibble(
    x = example_genotypes,
    loci = example_loci,
    indiv_meta = example_indiv_meta,
    quiet = TRUE
  )

  # gen_tibble with no missingness runs
  test_pca <- test %>% gt_pca_autoSVD(verbose = FALSE)

  # Now try with imputed data
  bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
  missing_gt <- gen_tibble(
    bed_file,
    backingfile = tempfile("missing_"),
    quiet = TRUE
  )
  expect_error(
    missing_gt %>% gt_pca_autoSVD(),
    "You can't have missing"
  )
  missing_gt <- gt_impute_simple(missing_gt, method = "mode")

  # we encounter our roll_size error
  expect_error(
    missing_pca <- missing_gt %>% gt_pca_autoSVD(verbose = FALSE),
    "roll_size exceeds the number of variants"
  )

  # adjusting roll_size fixes the error
  test_pca_roll0 <- missing_gt %>%
    gt_pca_autoSVD(roll_size = 0, verbose = FALSE)
  expect_s3_class(test_pca_roll0, "gt_pca")
  test_pca_roll7 <- missing_gt %>%
    gt_pca_autoSVD(roll_size = 7, verbose = FALSE)
  expect_s3_class(test_pca_roll7, "gt_pca")

  # testing with the families dataset too
  bed_file <- system.file(
    "extdata/related",
    "families.bed",
    package = "tidypopgen"
  )
  families <- gen_tibble(
    bed_file,
    backingfile = tempfile("families"),
    quiet = TRUE,
    valid_alleles = c("2", "1")
  )
  expect_error(
    families %>% gt_pca_autoSVD(),
    "You can't have missing"
  )
  families <- gt_impute_simple(families, method = "mode")

  # the same error occurs
  expect_error(
    missing_pca <- families %>% gt_pca_autoSVD(verbose = FALSE),
    "Parameter 'size' is too large."
  )

  # adjusting roll_size fixes the error
  test_pca_families_roll7 <- families %>%
    gt_pca_autoSVD(roll_size = 7, verbose = FALSE)
  expect_s3_class(test_pca_families_roll7, "gt_pca")
})


test_that("fit_gt_pca_and_predict_splitted_data", {
  bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
  missing_gt <- gen_tibble(
    bed_file,
    backingfile = tempfile("missing_"),
    quiet = TRUE
  )
  # create a fake ancient set by subsetting
  ancient_gt <- missing_gt[1:20, ]
  # now extract the modern data (to be imputed)
  modern_gt <- missing_gt[-c(1:20), ]

  # update the backingfiles
  ancient_gt <- gt_update_backingfile(ancient_gt, quiet = TRUE)
  modern_gt <- gt_update_backingfile(modern_gt, quiet = TRUE)

  # impute the modern data
  modern_gt <- gt_impute_simple(modern_gt, method = "mode")
  modern_pca <- modern_gt %>% gt_pca_partialSVD()
  # if we just try to predict, we find that the new data have missing data
  new_pred <- predict(
    modern_pca,
    new_data = ancient_gt,
    project_method = "simple",
    as_matrix = TRUE
  )
  expect_true(all(dim(new_pred) == c(20, 10)))
  # now raise an error if we don't impute to the mean
  expect_error(
    predict(modern_pca, new_data = ancient_gt, project_method = "none"),
    "You can't have missing values in 'X'"
  )
  # check the same with a tibble
  new_pred_tbl <- predict(
    modern_pca,
    new_data = ancient_gt,
    project_method = "simple",
    as_matrix = FALSE
  )
  expect_true(all(dim(new_pred_tbl) == c(20, 11)))
  expect_true(inherits(new_pred_tbl, "tbl_df"))
  # least squares prediction
  lsq_pred <- predict(
    modern_pca,
    new_data = ancient_gt,
    project_method = "least_squares",
    as_matrix = TRUE
  )
  expect_true(all(dim(lsq_pred) == c(20, 2)))
  lsq_pred_tbl <- predict(
    modern_pca,
    new_data = ancient_gt,
    project_method = "least_squares",
    as_matrix = FALSE
  )
  expect_true(all(dim(lsq_pred_tbl) == c(20, 3)))
  expect_true(inherits(lsq_pred_tbl, "tbl_df"))

  # test 3 components
  multi_comp_pred <- predict(
    modern_pca,
    new_data = ancient_gt,
    project_method = "least_squares",
    lsq_pcs = c(1, 2, 3),
    as_matrix = TRUE
  )
  expect_true(all(dim(multi_comp_pred) == c(20, 3)))

  # test errors from incorrect lsq_pcs
  expect_error(
    predict(
      modern_pca,
      new_data = ancient_gt,
      project_method = "least_squares",
      lsq_pcs = c(0, 2)
    ), "lsq_pcs should be a vector of valid"
  )
  expect_error(
    predict(
      modern_pca,
      new_data = ancient_gt,
      project_method = "least_squares",
      lsq_pcs = c(1, 20)
    ), "lsq_pcs should be a vector of valid"
  )

  # Test empty vector
  expect_error(
    predict(
      modern_pca,
      new_data = ancient_gt,
      project_method = "least_squares",
      lsq_pcs = c()
    ), "lsq_pcs should be a vector of valid"
  )

  # Test duplicate components
  expect_error(multi_comp_dup <- predict(
    modern_pca,
    new_data = ancient_gt,
    project_method = "least_squares",
    lsq_pcs = c(1, 1, 2),
    as_matrix = TRUE
  ), "lsq_pcs should not contain duplicate values")

  # Test single component
  single_comp <- predict(
    modern_pca,
    new_data = ancient_gt,
    project_method = "least_squares",
    lsq_pcs = 1,
    as_matrix = TRUE
  )
  expect_true(all(dim(single_comp) == c(20, 1)))

  # Test all components
  all_comp <- predict(
    modern_pca,
    new_data = ancient_gt,
    project_method = "least_squares",
    lsq_pcs = 1:10,
    as_matrix = TRUE
  )
  expect_true(all(dim(all_comp) == c(20, 10)))
})

test_that("PCA functions work with loci out of order", {
  bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
  missing_gt <- gen_tibble(
    bed_file,
    backingfile = tempfile("missing_"),
    quiet = TRUE
  )
  missing_gt <- gt_impute_simple(missing_gt, method = "mode")
  missing_part_pca1 <- missing_gt %>% gt_pca_partialSVD()
  missing_rand_pca1 <- missing_gt %>% gt_pca_randomSVD(verbose = FALSE)

  # now shuffle the loci
  loci <- missing_gt %>% show_loci()
  loci <- loci[sample(nrow(loci)), ]
  show_loci(missing_gt) <- loci

  # Rerun PCA
  missing_part_pca2 <- missing_gt %>% gt_pca_partialSVD()
  missing_rand_pca2 <- missing_gt %>% gt_pca_randomSVD(verbose = FALSE)

  expect_equal(missing_part_pca1$loadings, missing_part_pca2$loadings)
  expect_equal(missing_rand_pca1$loadings, missing_rand_pca2$loadings)
})

test_that("PCA computes frobenius when needed", {
  bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
  missing_gt <- gen_tibble(
    bed_file,
    backingfile = tempfile("missing_"),
    quiet = TRUE
  )
  missing_gt <- gt_impute_simple(missing_gt, method = "mode")
  missing_part_pca1 <- missing_gt %>% gt_pca_partialSVD()
  expect_true("square_frobenius" %in% names(missing_part_pca1))
  tidy_pca_out <- tidy(missing_part_pca1, matrix = "eigenvalues")
  expect_true("cumulative" %in% names(tidy_pca_out))
  # now repeat without estimating variance
  missing_part_pca2 <- missing_gt %>% gt_pca_partialSVD(total_var = FALSE)
  expect_false("square_frobenius" %in% names(missing_part_pca2))
  tidy_pca_out <- tidy(missing_part_pca2, matrix = "eigenvalues")
  expect_false("cumulative" %in% names(tidy_pca_out))
  # TODO we should repeat the PCA on a toy dataset to check that it is
  # equivalent to a full PCA done by hand
})

test_that("our stdevs are comparable to prcomp", {
  bed_path <- system.file(
    "extdata/related/families.bed",
    package = "tidypopgen"
  )
  families_bigsnp_path <- bigsnpr::snp_readBed(
    bed_path,
    backingfile = tempfile()
  )
  families <- gen_tibble(
    families_bigsnp_path,
    quiet = TRUE,
    valid_alleles = c("1", "2")
  )
  count_loci(families) # 941
  families <- families %>% select_loci_if(loci_maf(genotypes) > 0.01)
  # remove NA values for prcomp
  cols_with_na <- apply(show_genotypes(families), 2, function(x) any(is.na(x)))
  cols_without_na <- which(cols_with_na == FALSE)
  families <- families %>% select_loci(all_of(cols_without_na))
  count_loci(families) # 609

  #########################
  # Perform PCA using gt_pca_partialSVD
  #########################
  # TODO there seems to be a bug in the gt_pca_partialSVD function where
  # missing values are detected even when excluded
  families <- gt_impute_simple(families, method = "mode")
  gt_pca_result <- families %>% gt_pca_partialSVD()
  tidy_pca <- tidy(gt_pca_result, matrix = "eigenvalues")

  # Perform PCA using prcomp
  pca_result <- prcomp(
    show_genotypes(families),
    center = gt_pca_result$center,
    scale. = gt_pca_result$scale
  )
  prcomp_summary <- summary(pca_result)
  TOL <- 1e-4 # nolint

  # Compare standard deviation
  expect_equal(tidy_pca$std.dev, pca_result$sdev[1:10], tolerance = TOL)
  # Compare percentage variance explained
  expect_equal(
    tidy_pca$percent,
    as.vector((prcomp_summary$importance[2, c(1:10)]) * 100),
    tolerance = TOL
  )
  # Compare cumulative variance explained
  expect_equal(
    tidy_pca$cumulative,
    as.vector((prcomp_summary$importance[3, c(1:10)]) * 100),
    tolerance = TOL
  )

  #########################
  # Perform PCA using gt_pca_autoSVD
  #########################
  gt_pca_auto_result <- families %>% gt_pca_autoSVD(
    roll_size = 7,
    verbose = FALSE
  )
  tidy_pca <- tidy(gt_pca_auto_result, matrix = "eigenvalues")

  loci <- gt_pca_auto_result$loci
  select <- which(show_loci(families)$name %in% loci$name)
  families_autoSVD_subset <- families %>% # nolint
    select_loci(all_of(select)) %>%
    show_genotypes()

  # Perform PCA using prcomp
  pca_result <- prcomp(
    families_autoSVD_subset,
    center = gt_pca_auto_result$center,
    scale. = gt_pca_auto_result$scale
  )
  prcomp_summary <- summary(pca_result)

  # Compare standard deviation
  expect_equal(tidy_pca$std.dev, pca_result$sdev[1:10], tolerance = TOL)
  # Compare percentage variance explained
  expect_equal(
    tidy_pca$percent,
    as.vector((prcomp_summary$importance[2, c(1:10)]) * 100),
    tolerance = TOL
  )
  # Compare cumulative variance explained
  expect_equal(
    tidy_pca$cumulative,
    as.vector((prcomp_summary$importance[3, c(1:10)]) * 100),
    tolerance = TOL
  )
})

# The tests below require more than 2 cores, so we skip on CRAN
skip_on_cran()

test_that("n_cores can be set for pca functions", {
  skip_if(
    Sys.getenv("_R_CHECK_LIMIT_CORES_", "") != "",
    "Skipping due to core limitation in R CMD check"
  )
  bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
  missing_gt <- gen_tibble(
    bed_file,
    backingfile = tempfile("missing_"),
    quiet = TRUE
  )
  missing_gt <- gt_impute_simple(missing_gt, method = "mode")

  ############
  # Test gt_pca_randomSVD
  ############
  expect_true(getOption("bigstatsr.check.parallel.blas"))
  one_core <- missing_gt %>%
    gt_pca_randomSVD(n_cores = 1, verbose = FALSE)
  two_core <- missing_gt %>%
    gt_pca_randomSVD(n_cores = 2, verbose = FALSE)
  # calls will differ, set to NULL
  one_core$call <- NULL
  two_core$call <- NULL
  expect_equal(one_core, two_core)
  expect_true(getOption("bigstatsr.check.parallel.blas"))

  # test parallel blas is true on exit if function errors
  expect_error(
    gt_pca_autoSVD(one_core, n_cores = 2),
    "no applicable method for 'show_loci'"
  )
  expect_true(getOption("bigstatsr.check.parallel.blas"))
})

test_that("n_cores can be set gt_pca_autoSVD", {
  skip_if(
    Sys.getenv("_R_CHECK_LIMIT_CORES_", "") != "",
    "Skipping due to core limitation in R CMD check"
  )
  bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
  missing_gt <- gen_tibble(
    bed_file,
    backingfile = tempfile("missing_"),
    quiet = TRUE
  )
  missing_gt <- gt_impute_simple(missing_gt, method = "mode")
  ############
  # Test gt_pca_autoSVD
  ############
  expect_true(getOption("bigstatsr.check.parallel.blas"))
  one_core <- missing_gt %>%
    gt_pca_autoSVD(n_cores = 1, roll_size = 7, verbose = FALSE)
  two_core <- missing_gt %>%
    gt_pca_autoSVD(n_cores = 2, roll_size = 7, verbose = FALSE)
  # calls will differ, set to NULL
  one_core$call <- NULL
  two_core$call <- NULL
  expect_equal(one_core, two_core)
  expect_true(getOption("bigstatsr.check.parallel.blas"))

  # test parallel blas is true on exit if function errors
  expect_error(
    gt_pca_autoSVD(one_core, n_cores = 2),
    "no applicable method for 'show_loci'"
  )
  expect_true(getOption("bigstatsr.check.parallel.blas"))
})


test_that("n_cores can be set for predict_gt_pca", {
  bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
  missing_gt <- gen_tibble(
    bed_file,
    backingfile = tempfile("missing_"),
    quiet = TRUE
  )
  # create a fake ancient set by subsetting
  ancient_gt <- missing_gt[1:20, ]
  # now extract the modern data (to be imputed)
  modern_gt <- missing_gt[-c(1:20), ]

  # update the backingfiles
  ancient_gt <- gt_update_backingfile(ancient_gt, quiet = TRUE)
  modern_gt <- gt_update_backingfile(modern_gt, quiet = TRUE)

  modern_gt <- gt_impute_simple(modern_gt, method = "mode")
  modern_pca <- modern_gt %>% gt_pca_partialSVD()
  expect_true(getOption("bigstatsr.check.parallel.blas"))
  one_core <- predict(
    modern_pca,
    new_data = ancient_gt,
    project_method = "simple",
    n_cores = 1
  )
  two_core <- predict(
    modern_pca,
    new_data = ancient_gt,
    project_method = "simple",
    n_cores = 2
  )
  expect_equal(one_core, two_core)
  expect_true(getOption("bigstatsr.check.parallel.blas"))

  # test parallel blas is true on exit if function errors
  expect_error(
    predict(modern_pca, new_data = ancient_gt, project_method = "none"),
    "You can't have missing values in 'X'"
  )
  expect_true(getOption("bigstatsr.check.parallel.blas"))
})
