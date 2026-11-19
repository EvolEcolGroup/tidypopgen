options(mc_doScale_quiet = TRUE)

test_that("fit_gt_pca_and_predict", {
  # the deprecated gt_pca_partialSVD() errors on missing data
  withr::local_options(lifecycle_verbosity = "quiet")
  bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
  missing_gt <- gen_tibble(
    bed_file,
    backingfile = tempfile("missing_"),
    quiet = TRUE
  )
  expect_error(
    missing_gt %>% gt_pca_partialSVD(),
    "needs complete or imputed genotypes"
  )
  missing_gt <- gt_impute_simple(missing_gt, method = "mode")
  # the pca of the imputed data, which predict() uses for new_data, in double
  # precision (the internal engine's single precision products differ from
  # the projection by about 1e-7)
  gt_set_imputed(missing_gt, set = TRUE)
  missing_pca <- missing_gt %>% gt_pca(n_cores = 2, engine = "blas")
  gt_set_imputed(missing_gt, set = FALSE)
  # check that predicting on the object is the same as predicting from the
  # full dataset without imputation to the center (the data are already imputed)
  # ("none" uses the imputed genotypes automatically; the default least
  # squares only if the gen_tibble is set to use them)
  expect_true(all.equal(
    predict(missing_pca),
    predict(missing_pca, new_data = missing_gt, project_method = "none"),
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
  gt_set_imputed(missing_gt_sub, set = TRUE)
  missing_sub_pca <- missing_gt_sub %>% gt_pca(n_cores = 2, engine = "blas")
  gt_set_imputed(missing_gt_sub, set = FALSE)
  expect_true(all.equal(
    predict(missing_sub_pca),
    predict(missing_sub_pca, new_data = missing_gt, project_method = "none"),
    check.attributes = FALSE
  ))
})

test_that("adjusting roll_size fixes gt_pca_autoSVD problem ", {
  # the genotypes are random: without a seed, the first gt_pca_autoSVD()
  # call occasionally failed with the roll_size error
  set.seed(1)
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
  # and records complete genotypes, not imputed
  expect_false(attr(test_pca, "imputed"))
  expect_false(attr(test_pca, "missing"))

  # Now try with imputed data
  bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
  missing_gt <- gen_tibble(
    bed_file,
    backingfile = tempfile("missing_"),
    quiet = TRUE
  )
  expect_error(
    missing_gt %>% gt_pca_autoSVD(),
    "needs complete or imputed genotypes"
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
  # the imputed genotypes were used
  expect_true(attr(test_pca_roll7, "imputed"))
  expect_false(attr(test_pca_roll7, "missing"))

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
    "needs complete or imputed genotypes"
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
  modern_pca <- modern_gt %>% gt_pca(n_cores = 2)
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
    "needs complete or imputed genotypes"
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
  # all the components by default
  expect_true(all(dim(lsq_pred) == c(20, ncol(modern_pca$u))))
  lsq_pred_tbl <- predict(
    modern_pca,
    new_data = ancient_gt,
    project_method = "least_squares",
    as_matrix = FALSE
  )
  expect_true(all(dim(lsq_pred_tbl) == c(20, ncol(modern_pca$u) + 1)))
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

  # Test empty vector (NULL selects all the components)
  expect_error(
    predict(
      modern_pca,
      new_data = ancient_gt,
      project_method = "least_squares",
      lsq_pcs = integer(0)
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
  missing_part_pca1 <- missing_gt %>% gt_pca(n_cores = 2)
  missing_rand_pca1 <- missing_gt %>%
    gt_pca(method = "random", n_cores = 1, verbose = FALSE)

  # now shuffle the loci
  loci <- missing_gt %>% show_loci()
  loci <- loci[sample(nrow(loci)), ]
  show_loci(missing_gt) <- loci

  # Rerun PCA
  missing_part_pca2 <- missing_gt %>% gt_pca(n_cores = 2)
  missing_rand_pca2 <- missing_gt %>%
    gt_pca(method = "random", n_cores = 1, verbose = FALSE)

  expect_equal(missing_part_pca1$loadings, missing_part_pca2$loadings)
  expect_equal(missing_rand_pca1$loadings, missing_rand_pca2$loadings)
})

test_that("gt_pca_autoSVD works with loci out of order", {
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
  pca_ordered <- gt_pca_autoSVD(test_gt,
    k = 3, use_positions = FALSE,
    roll_size = 7, verbose = FALSE
  )
  # reverse the order of the loci in the loci table
  test_gt_rev <- test_gt %>% select_loci(rev(seq_len(n_loci)))
  pca_rev <- gt_pca_autoSVD(test_gt_rev,
    k = 3, use_positions = FALSE,
    roll_size = 7, verbose = FALSE
  )
  # the same loci should be retained, and the loci table should be in the
  # same order as the loadings
  expect_identical(pca_rev$loci$name, row_names(pca_rev$v))
  expect_identical(pca_ordered$loci$name, row_names(pca_ordered$v))
  expect_identical(sort(pca_rev$loci$name), sort(pca_ordered$loci$name))
  # loadings, centre and scale of each locus are the same
  match_loci <- match(row_names(pca_ordered$v), row_names(pca_rev$v))
  expect_equal(
    abs(unclass(pca_rev$v)[match_loci, ]),
    abs(unclass(pca_ordered$v)),
    ignore_attr = TRUE
  )
  expect_equal(pca_rev$center[match_loci], pca_ordered$center)
  expect_equal(pca_rev$scale[match_loci], pca_ordered$scale)
  # the total variance is the same
  expect_equal(pca_rev$square_frobenius, pca_ordered$square_frobenius)
  expect_equal(
    tidy(pca_rev, "eigenvalues")$percent,
    tidy(pca_ordered, "eigenvalues")$percent
  )
  # predicting the original data recovers the scores
  expect_equal(
    abs(predict(pca_rev, new_data = test_gt_rev, project_method = "none")),
    abs(predict(pca_rev)),
    ignore_attr = TRUE
  )
  expect_equal(
    abs(predict(pca_rev, new_data = test_gt, project_method = "none")),
    abs(predict(pca_ordered)),
    ignore_attr = TRUE
  )
})

test_that("PCA computes frobenius when needed", {
  bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
  missing_gt <- gen_tibble(
    bed_file,
    backingfile = tempfile("missing_"),
    quiet = TRUE
  )
  missing_gt <- gt_impute_simple(missing_gt, method = "mode")
  missing_part_pca1 <- missing_gt %>% gt_pca(n_cores = 2)
  expect_true("square_frobenius" %in% names(missing_part_pca1))
  tidy_pca_out <- tidy(missing_part_pca1, matrix = "eigenvalues")
  expect_true("cumulative" %in% names(tidy_pca_out))
  # now repeat without estimating variance
  missing_part_pca2 <- missing_gt %>% gt_pca(total_var = FALSE, n_cores = 2)
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
  # Perform PCA using gt_pca
  #########################
  families <- gt_impute_simple(families, method = "mode")
  gt_pca_result <- families %>% gt_pca(n_cores = 2)
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
  # Test gt_pca with method "random"
  ############
  expect_true(getOption("bigstatsr.check.parallel.blas"))
  one_core <- missing_gt %>%
    gt_pca(method = "random", n_cores = 1, verbose = FALSE)
  two_core <- missing_gt %>%
    gt_pca(method = "random", n_cores = 2, verbose = FALSE)
  # calls will differ, set to NULL
  one_core$call <- NULL
  two_core$call <- NULL
  expect_equal(one_core, two_core)
  expect_true(getOption("bigstatsr.check.parallel.blas"))

  # test parallel blas is true on exit if function errors (and bigstatsr's
  # checks are turned on again)
  expect_error(
    gt_pca_autoSVD(one_core, n_cores = 2),
    "not a gen_tibble"
  )
  expect_true(getOption("bigstatsr.check.parallel.blas"))
  expect_true(getOption("bigstatsr.check.args"))
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

  # test parallel blas is true on exit if function errors (and bigstatsr's
  # checks are turned on again)
  expect_error(
    gt_pca_autoSVD(one_core, n_cores = 2),
    "not a gen_tibble"
  )
  expect_true(getOption("bigstatsr.check.parallel.blas"))
  expect_true(getOption("bigstatsr.check.args"))
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
  modern_pca <- modern_gt %>% gt_pca(n_cores = 2)
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
    "needs complete or imputed genotypes"
  )
  expect_true(getOption("bigstatsr.check.parallel.blas"))
})

test_that("gt_pca matches PLINK 1.9 with very high missingness", {
  # 60 individuals with 0 to 90% missing genotypes, and the outputs of
  # plink --make-rel square and plink --pca 5
  # (made by data-raw/validate/create_plink_pca_test_data.R)
  fixture <- test_path("testdata", "plink_pca", "missing")
  missing_gt <- gen_tibble(paste0(fixture, ".bed"),
    backingfile = tempfile("grm_"), quiet = TRUE
  )
  # the relationship matrix, against PLINK and an R reference
  freq <- loci_alt_freq(missing_gt)
  geno <- show_genotypes(missing_gt)
  z <- sweep(sweep(geno, 2, 2 * freq), 2, sqrt(2 * freq * (1 - freq)), "/")
  z[!is.finite(z)] <- 0 # missing genotypes and the monomorphic locus
  valid <- (!is.na(geno)) * 1
  grm_ref <- tcrossprod(z) / tcrossprod(valid)
  # both engines of the matrix products: the BLAS (double precision), and
  # Eigen (single precision products of up to 1,024 loci, summed in double)
  for (use_blas in c(TRUE, FALSE)) {
    grm <- grm_plink_cpp(
      BM = .gt_get_fbm(missing_gt), rowInd = .gt_fbm_rows(missing_gt),
      colInd = show_loci(missing_gt)$big_index, center = 2 * freq,
      scale = sqrt(2 * freq * (1 - freq)), block_size = 64, ncores = 2,
      use_blas = use_blas
    )
    # the number of missing genotypes (for the attributes of gt_pca objects)
    expect_equal(attr(grm, "n_missing"), sum(is.na(geno)))
    attr(grm, "n_missing") <- NULL
    expect_equal(grm, grm_ref, tolerance = if (use_blas) 1e-12 else 1e-6)
  }
  plink_rel <- unname(as.matrix(read.table(paste0(fixture, ".rel"))))
  # PLINK writes 6 significant digits
  expect_lt(max(abs(grm - plink_rel)), 1e-5)
  # the principal components, against PLINK (up to their sign)
  missing_pca <- gt_pca(missing_gt, k = 5, n_cores = 2)
  plink_val <- scan(paste0(fixture, ".eigenval"), quiet = TRUE)
  plink_vec <- as.matrix(read.table(paste0(fixture, ".eigenvec"))[, -(1:2)])
  expect_equal(missing_pca$d^2 / count_loci(missing_gt), plink_val,
    tolerance = 1e-5
  )
  signs <- sign(colSums(missing_pca$u * plink_vec))
  expect_lt(max(abs(sweep(missing_pca$u, 2, signs, "*") - plink_vec)), 1e-5)
  expect_equal(rownames(missing_pca$u), missing_gt$id)
  # the result works with the gt_pca methods
  expect_s3_class(missing_pca, "gt_pca")
  expect_equal(nrow(tidy(missing_pca, "eigenvalues")), 5)
  expect_equal(dim(predict(missing_pca)), c(nrow(missing_gt), 5))
  # the same with reordered individuals and loci
  sub_gt <- missing_gt %>%
    dplyr::slice(60:1) %>%
    select_loci(1000:1)
  sub_pca <- gt_pca(sub_gt,
    k = 5, n_cores = 1, block_size = 7, engine = "internal"
  )
  # (the internal engine's single precision products differ by about 1e-8)
  expect_equal(sub_pca$d, missing_pca$d, tolerance = 1e-6)
  expect_equal(abs(sub_pca$u[60:1, ]), abs(missing_pca$u),
    tolerance = 1e-6, ignore_attr = TRUE
  )
})

test_that("predict(project_method = 'plink') matches PLINK 1.9", {
  # PLINK's PCA of the 24 individuals with 0 or 20% missing genotypes, and its
  # projection of the other 36 (50 to 90% missing), with --pca-cluster-names
  # (made by data-raw/validate/create_plink_pca_test_data.R)
  fixture <- test_path("testdata", "plink_pca", "missing")
  missing_gt <- gen_tibble(paste0(fixture, ".bed"),
    backingfile = tempfile("proj_"), quiet = TRUE
  )
  # PLINK drops the monomorphic loci with --pca-cluster-names
  missing_gt <- select_loci_if(missing_gt, loci_maf(genotypes) > 0)
  is_ref <- rep(rep(c(TRUE, FALSE), c(8, 12)), 3)
  # PLINK standardises with the frequencies of all the individuals
  freq <- loci_alt_freq(missing_gt)
  all_freq <- function(...) {
    data.frame(center = 2 * freq, scale = sqrt(2 * freq * (1 - freq)))
  }
  ref_pca <- gt_pca(missing_gt[is_ref, ],
    k = 3, fun_scaling = all_freq, n_cores = 1
  )
  plink_vec <- as.matrix(
    read.table(paste0(fixture, "_proj.eigenvec"))[, -(1:2)]
  )
  signs <- sign(colSums(ref_pca$u * plink_vec[is_ref, ]))
  expect_lt(
    max(abs(sweep(ref_pca$u, 2, signs, "*") - plink_vec[is_ref, ])), 1e-5
  )
  # PLINK's eigenvectors of projected individuals are scores / d
  plink_proj <- predict(ref_pca,
    new_data = missing_gt[!is_ref, ], project_method = "plink"
  )
  plink_scaled <- sweep(plink_proj, 2, signs * ref_pca$d, "/")
  expect_lt(max(abs(plink_scaled - plink_vec[!is_ref, ])), 1e-5)
})

test_that("gt_pca matches gt_pca_partialSVD on complete data", {
  # the deprecated function is the reference
  withr::local_options(lifecycle_verbosity = "quiet")
  bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
  complete_gt <- gen_tibble(bed_file,
    backingfile = tempfile("grm_complete_"), quiet = TRUE
  )
  complete_gt <- gt_impute_simple(complete_gt, method = "mode")
  gt_set_imputed(complete_gt, set = TRUE)
  # drop monomorphic loci, which partialSVD can not scale
  freq <- loci_alt_freq(complete_gt)
  complete_gt <- select_loci(complete_gt, which(freq > 0 & freq < 1))
  svd_pca <- gt_pca_partialSVD(complete_gt, k = 4)
  # both engines explicitly: "auto" picks the internal engine with R's
  # reference BLAS (Windows, macOS) or when it times faster, and its single
  # precision products differ by about 1e-7
  for (engine in c("blas", "internal")) {
    tol <- if (engine == "blas") 1e-8 else 1e-6
    grm_pca <- gt_pca(complete_gt, k = 4, engine = engine)
    expect_equal(grm_pca$d, svd_pca$d, tolerance = tol)
    expect_equal(abs(grm_pca$u), abs(svd_pca$u), tolerance = tol,
      ignore_attr = TRUE
    )
    expect_equal(abs(grm_pca$v), abs(svd_pca$v), tolerance = tol,
      ignore_attr = TRUE
    )
  }
  gt_set_imputed(complete_gt, set = FALSE)
})

test_that("gt_pca rejects pairs without shared loci and bad k", {
  geno <- rbind(c(0, 1, NA, NA), c(NA, NA, 2, 1), c(1, 2, 0, 1))
  test_gt <- gen_tibble(geno,
    indiv_meta = data.frame(id = c("a", "b", "c"), population = "p"),
    loci = data.frame(
      name = paste0("rs", 1:4), chromosome = 1, position = 1:4,
      genetic_dist = 0, allele_ref = "A", allele_alt = "T"
    ),
    backingfile = tempfile("grm_small_"), quiet = TRUE
  )
  expect_error(gt_pca(test_gt, k = 1), "no locus genotyped in both")
  expect_error(gt_pca(test_gt, k = 3), "'k' should be between")
})

test_that("gt_pca method random needs complete or imputed genotypes", {
  withr::local_options(lifecycle_verbosity = "quiet")
  set.seed(1)
  sim_gt <- gt_simulate_bn(
    n_indiv = c(20, 20), n_loci = 500, fst = 0.1, missing = 0.1
  )
  expect_error(gt_pca(sim_gt, k = 2, method = "random"), "needs complete")
  # with imputed genotypes, used automatically, as gt_pca_randomSVD() does
  sim_gt <- gt_impute_simple(sim_gt, method = "mode")
  gt_set_imputed(sim_gt, set = TRUE)
  sim_gt <- select_loci_if(sim_gt, loci_maf(genotypes) > 0)
  gt_set_imputed(sim_gt, set = FALSE)
  random_pca <- gt_pca(sim_gt, k = 2, method = "random", n_cores = 1)
  old_pca <- gt_pca_randomSVD(sim_gt, k = 2)
  expect_equal(random_pca$d, old_pca$d)
  expect_equal(random_pca$method, "random")
  expect_false(gt_uses_imputed(sim_gt))
  # and the same as the partial method on the imputed genotypes
  gt_set_imputed(sim_gt, set = TRUE)
  partial_pca <- gt_pca(sim_gt, k = 2, n_cores = 1)
  gt_set_imputed(sim_gt, set = FALSE)
  expect_equal(random_pca$d, partial_pca$d, tolerance = 1e-4)
})

test_that("missing genotypes are looked for only where they are used", {
  withr::local_options(lifecycle_verbosity = "quiet")
  # example-missing has missing genotypes at 367 of its 500 loci: keep the
  # 133 complete loci, without resaving the backing matrix (where the other
  # loci keep their missing genotypes)
  bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
  missing_gt <- gen_tibble(
    bed_file,
    backingfile = tempfile("missing_"),
    quiet = TRUE
  )
  complete_gt <- missing_gt %>%
    select_loci_if(loci_missingness(genotypes) == 0)
  expect_false(any(loci_missingness(complete_gt) > 0))
  # bigstatsr's check sampled the whole backing matrix, and stopped with
  # "You can't have missing values" (for the workers of the parallel
  # randomised SVD, see the next test)
  set.seed(1)
  expect_s3_class(gt_pca_partialSVD(complete_gt, k = 2), "gt_pca")
  expect_s3_class(
    gt_pca(complete_gt, k = 2, method = "random", n_cores = 1), "gt_pca"
  )
  expect_s3_class(
    gt_pca_autoSVD(complete_gt,
      k = 2, roll_size = 0, verbose = FALSE, n_cores = 1
    ),
    "gt_pca"
  )
  complete_pca <- gt_pca(complete_gt, k = 2, n_cores = 1)
  expect_identical(
    dim(predict(complete_pca, new_data = complete_gt, project_method = "none")),
    c(nrow(complete_gt), 2L)
  )
  # missing genotypes in the loci used are still found
  expect_error(
    gt_pca_partialSVD(missing_gt, k = 2),
    "gt_pca_partialSVD\\(\\) needs complete or imputed genotypes"
  )
  expect_error(
    gt_pca_randomSVD(missing_gt, k = 2),
    "gt_pca_randomSVD\\(\\) needs complete or imputed genotypes"
  )
})

test_that("missing genotypes not used are ignored by the parallel SVD", {
  # the parallel randomised SVD runs 1 + n_cores processes, more than the 2
  # allowed by R CMD check
  skip_if(
    Sys.getenv("_R_CHECK_LIMIT_CORES_", "") != "",
    "Skipping due to core limitation in R CMD check"
  )
  withr::local_options(lifecycle_verbosity = "quiet")
  bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
  missing_gt <- gen_tibble(
    bed_file,
    backingfile = tempfile("missing_"),
    quiet = TRUE
  )
  complete_gt <- missing_gt %>%
    select_loci_if(loci_missingness(genotypes) == 0)
  # bigstatsr's check in the workers sampled the whole backing matrix, and
  # left the SVD waiting forever
  set.seed(1)
  expect_s3_class(
    gt_pca(complete_gt, k = 2, method = "random", n_cores = 2), "gt_pca"
  )
  expect_s3_class(
    gt_pca_autoSVD(complete_gt,
      k = 2, roll_size = 0, verbose = FALSE, n_cores = 2
    ),
    "gt_pca"
  )
})

test_that("the SVDs check all the genotypes they use, also in parallel", {
  withr::local_options(lifecycle_verbosity = "quiet")
  set.seed(6)
  sim_gt <- gt_simulate_bn(n_indiv = c(10, 10), n_loci = 300, fst = 0.1)
  sim_gt <- select_loci_if(sim_gt, loci_maf(genotypes) > 0)
  # a single missing genotype is found by the scaling (computed once, before
  # the SVD), which gave a cryptic error, or, with n_cores > 1, left the
  # workers of the randomised SVD waiting forever
  one_missing <- show_genotypes(sim_gt)
  one_missing[3, 150] <- NA
  one_missing_gt <- gen_tibble(one_missing,
    indiv_meta = sim_gt %>% select(id, population),
    loci = show_loci(sim_gt) %>% select(-big_index),
    backingfile = tempfile(), quiet = TRUE
  )
  expect_error(
    gt_pca(one_missing_gt, k = 2, method = "random", n_cores = 2),
    "gt_pca\\(method = \"random\"\\) needs complete or imputed genotypes"
  )
  expect_error(
    gt_pca_autoSVD(one_missing_gt,
      k = 2, roll_size = 0, verbose = FALSE, n_cores = 2
    ),
    "gt_pca_autoSVD\\(\\) needs complete or imputed genotypes"
  )
  # with a scaling that does not depend on the genotypes, the sampled check
  # finds it (the SVD failed with "TridiagEigen: eigen decomposition failed")
  constant_scaling <- bigstatsr::big_scale(center = FALSE, scale = FALSE)
  expect_error(
    gt_pca(one_missing_gt,
      k = 2, method = "random", fun_scaling = constant_scaling, n_cores = 1
    ),
    "gt_pca\\(method = \"random\"\\) needs complete or imputed genotypes"
  )
  expect_error(
    gt_pca_partialSVD(one_missing_gt, k = 2, fun_scaling = constant_scaling),
    "gt_pca_partialSVD\\(\\) needs complete or imputed genotypes"
  )
  # and monomorphic loci give bigstatsr's error before the SVD starts
  mono_gt <- gt_simulate_bn(n_indiv = c(10, 10), n_loci = 300, fst = 0.1)
  expect_true(any(loci_maf(mono_gt) == 0))
  expect_error(
    gt_pca(mono_gt, k = 2, method = "random", n_cores = 2),
    "zero scaling"
  )
})

test_that("the number of cores is checked, as bigstatsr did", {
  set.seed(7)
  sim_gt <- gt_simulate_bn(n_indiv = c(10, 10), n_loci = 200, fst = 0.1)
  sim_gt <- select_loci_if(sim_gt, loci_maf(genotypes) > 0)
  # (bigstatsr's own checks are off, for the genotypes; this one is kept)
  expect_error(
    gt_pca(sim_gt, k = 2, method = "random", n_cores = 0),
    "'ncores' should be at least 1"
  )
  expect_error(
    gt_pca_autoSVD(sim_gt, k = 2, roll_size = 0, verbose = FALSE, n_cores = 0),
    "'ncores' should be at least 1"
  )
  expect_error(
    gt_pca(sim_gt,
      k = 2, method = "random",
      n_cores = getOption("bigstatsr.ncores.max") + 1
    ),
    "more cores than allowed"
  )
  expect_true(getOption("bigstatsr.check.args"))
})

test_that("gt_pca objects record imputed and missing genotypes", {
  withr::local_options(lifecycle_verbosity = "quiet")
  set.seed(1)
  sim_gt <- gt_simulate_bn(
    n_indiv = c(20, 20), n_loci = 500, fst = 0.1, missing = 0.1
  )
  sim_gt <- select_loci_if(sim_gt, loci_maf(genotypes) > 0)
  # partial, with missing genotypes compared pairwise
  missing_pca <- gt_pca(sim_gt, k = 2, n_cores = 1)
  expect_false(attr(missing_pca, "imputed"))
  expect_true(attr(missing_pca, "missing"))
  expect_output(print(missing_pca), "with missing genotypes")
  # partial, on imputed genotypes (only if the gen_tibble uses them)
  sim_gt <- gt_impute_simple(sim_gt, method = "mode")
  expect_true(attr(gt_pca(sim_gt, k = 2, n_cores = 1), "missing"))
  gt_set_imputed(sim_gt, set = TRUE)
  imputed_pca <- gt_pca(sim_gt, k = 2, n_cores = 1)
  gt_set_imputed(sim_gt, set = FALSE)
  expect_true(attr(imputed_pca, "imputed"))
  expect_false(attr(imputed_pca, "missing"))
  expect_output(print(imputed_pca), "Genotypes: imputed")
  # random and the deprecated functions use the imputed genotypes
  for (pca in list(
    gt_pca(sim_gt, k = 2, method = "random", n_cores = 1),
    gt_pca_randomSVD(sim_gt, k = 2),
    gt_pca_partialSVD(sim_gt, k = 2)
  )) {
    expect_true(attr(pca, "imputed"))
    expect_false(attr(pca, "missing"))
  }
  # partial, on complete genotypes
  complete_gt <- gt_simulate_bn(n_indiv = c(10, 10), n_loci = 200, fst = 0.1)
  complete_pca <- gt_pca(complete_gt, k = 2, n_cores = 1)
  expect_false(attr(complete_pca, "imputed"))
  expect_false(attr(complete_pca, "missing"))
  expect_output(print(complete_pca), "Genotypes: complete")
})

test_that("gt_pca accepts a custom scaling, as gt_pca_partialSVD does", {
  withr::local_options(lifecycle_verbosity = "quiet")
  set.seed(2)
  sim_gt <- gt_simulate_bn(n_indiv = c(15, 15), n_loci = 400, fst = 0.1)
  sim_gt <- select_loci_if(sim_gt, loci_maf(genotypes) > 0)
  for (scaling in list(bigstatsr::big_scale(), bigsnpr::snp_scaleAlpha(-0.5))) {
    new_pca <- gt_pca(sim_gt, k = 3, fun_scaling = scaling, n_cores = 1)
    old_pca <- gt_pca_partialSVD(sim_gt, k = 3, fun_scaling = scaling)
    expect_equal(new_pca$d, old_pca$d, tolerance = 1e-6)
    expect_equal(new_pca$center, old_pca$center)
    expect_equal(abs(new_pca$u), abs(old_pca$u),
      tolerance = 1e-6, ignore_attr = TRUE
    )
  }
  # a custom scaling can not be computed with missing genotypes
  sim_missing <- gt_simulate_bn(n_indiv = 20, n_loci = 300, missing = 0.2)
  expect_error(
    gt_pca(sim_missing, k = 2, fun_scaling = bigstatsr::big_scale()),
    "complete"
  )
})

test_that("the old PCA functions are deprecated", {
  withr::local_options(lifecycle_verbosity = "warning")
  set.seed(3)
  sim_gt <- gt_simulate_bn(n_indiv = c(10, 10), n_loci = 200, fst = 0.1)
  sim_gt <- select_loci_if(sim_gt, loci_maf(genotypes) > 0)
  expect_warning(gt_pca_partialSVD(sim_gt, k = 2), "deprecated")
  expect_warning(gt_pca_randomSVD(sim_gt, k = 2), "deprecated")
})

test_that("gt_pca works on a grouped gen_tibble", {
  set.seed(4)
  sim_gt <- gt_simulate_bn(n_indiv = c(10, 10), n_loci = 200, fst = 0.1)
  sim_gt <- select_loci_if(sim_gt, loci_maf(genotypes) > 0)
  grouped_gt <- sim_gt %>% group_by(population)
  # the grouping is ignored, with both methods
  expect_equal(
    gt_pca(grouped_gt, k = 2, n_cores = 1)$d,
    gt_pca(sim_gt, k = 2, n_cores = 1)$d
  )
  expect_equal(
    gt_pca(grouped_gt, k = 2, method = "random", n_cores = 1)$d,
    gt_pca(sim_gt, k = 2, method = "random", n_cores = 1)$d
  )
})

test_that("projections with missing genotypes match an R reference", {
  set.seed(5)
  sim_gt <- gt_simulate_bn(
    n_indiv = c(15, 15), n_loci = 300, fst = 0.1,
    missing = rep(c(0, 0.3, 0.6), 10)
  )
  sim_gt <- select_loci_if(sim_gt, loci_maf(genotypes) > 0)
  pca <- gt_pca(sim_gt[1:20, ], k = 3, n_cores = 1)
  new_gt <- sim_gt[21:30, ]
  # standardised genotypes of the new individuals, missing as 0 (and loci
  # monomorphic in the reference, with a scale of 0, also 0)
  geno <- show_genotypes(new_gt)
  z <- sweep(sweep(geno, 2, pca$center), 2, pca$scale, "/")
  z[!is.finite(z)] <- 0
  v <- unclass(pca$v)
  simple <- predict(pca,
    new_data = new_gt, project_method = "simple", block_size = 7
  )
  expect_equal(simple, z %*% v, ignore_attr = TRUE)
  # "none" on the complete individuals gives the same scores, also with loci
  # monomorphic in the reference (a scale of 0)
  expect_true(any(pca$scale == 0))
  # ("none" needs imputed genotypes, which change nothing for these)
  complete <- which(rowSums(is.na(geno)) == 0)
  complete_gt <- gt_update_backingfile(new_gt[complete, ], quiet = TRUE)
  complete_gt <- gt_impute_simple(complete_gt, method = "mode")
  none_proj <- predict(pca, new_data = complete_gt, project_method = "none")
  expect_equal(none_proj, simple[complete, ], ignore_attr = TRUE)
  # "plink": the simple projection divided by the fraction of loci genotyped
  plink_proj <- predict(pca,
    new_data = new_gt, project_method = "plink", block_size = 7
  )
  expect_equal(plink_proj, z %*% v / rowMeans(!is.na(geno)),
    ignore_attr = TRUE
  )
  # least squares on the loci genotyped in each individual
  lsq_ref <- t(sapply(seq_len(nrow(geno)), function(i) {
    obs <- !is.na(geno[i, ])
    solve(crossprod(v[obs, , drop = FALSE]), crossprod(v[obs, ], z[i, obs]))
  }))
  lsq <- predict(pca,
    new_data = new_gt, project_method = "least_squares", lsq_pcs = 1:3,
    block_size = 7
  )
  expect_equal(lsq, lsq_ref, ignore_attr = TRUE)
  # individuals genotyped at fewer loci than components get NA, with a warning
  set.seed(2)
  sparse_gt <- gt_simulate_bn(n_indiv = 6, n_loci = 300, missing = 0.99)
  # genotyped loci among those of the pca (predict() matches them by name)
  sparse_geno <- show_genotypes(sparse_gt)
  sparse_geno <- sparse_geno[, match(pca$loci$name, loci_names(sparse_gt))]
  n_obs <- rowSums(!is.na(sparse_geno))
  expect_warning(
    sparse <- predict(pca, new_data = sparse_gt, as_matrix = TRUE),
    "too few genotyped loci"
  )
  expect_true(any(n_obs < 3) && any(n_obs >= 3))
  expect_identical(is.na(sparse[, 1]), n_obs < 3, ignore_attr = TRUE)
})
