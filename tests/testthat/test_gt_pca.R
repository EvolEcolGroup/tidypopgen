options(mc_doScale_quiet=TRUE)

test_that("fit_gt_pca_and_predict",{
  bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
  missing_gt <- gen_tibble(bed_file,  backingfile = tempfile("missing_"),quiet=TRUE)
  expect_error( missing_gt %>% gt_pca_partialSVD(),
                "You can't have missing")
  missing_gt <- gt_impute_simple(missing_gt, method = "mode")
  missing_pca <- missing_gt %>% gt_pca_partialSVD()
  # check that predicting on the object is the same as predicting from the full dataset
  # without imputation to the center (the data are already imputed)
  expect_true(all.equal(predict(missing_pca),
                        predict(missing_pca, new_data = missing_gt),
                        check.attributes=FALSE))
  # now mismatch the loci table
  missing_gt_edited <- missing_gt
  show_loci(missing_gt_edited)$name[3] <- "blah"
  expect_error(predict(missing_pca, new_data = missing_gt_edited),
               "loci used in object")
  missing_gt_edited <- missing_gt
  show_loci(missing_gt_edited)$allele_ref[3] <- "blah"
  expect_error(predict(missing_pca, new_data = missing_gt_edited),
               "ref and alt alleles differ")
  # predict when new dataset has extra positions
  missing_gt_sub <- missing_gt %>% select_loci(100:450)
  missing_sub_pca <- missing_gt_sub  %>% gt_pca_partialSVD()
  expect_true(all.equal(predict(missing_sub_pca),
                        predict(missing_sub_pca, new_data = missing_gt,
                                project_method = "none"),
                        check.attributes=FALSE))

})

test_that("adjusting roll_size fixes gt_pca_autoSVD problem ",{

  # @TODO remove this hack when the patch makes it through to cran

  # Generate example_indiv_meta data frame
  individual_ids <- paste0("indiv", 1:200)
  populations <- sample(c("pop1", "pop2", "pop3", "pop4"), size = 200, replace = TRUE)
  example_indiv_meta <- data.frame(id = individual_ids, population = populations)

  # Generate example_genotypes matrix (no missingness)
  values <- c(0, 1, 2)
  genotype_matrix <- sample(values, size = 200 * 500, replace = TRUE)
  example_genotypes <- matrix(genotype_matrix, nrow = 200, ncol = 500)

  # Generate example_loci data frame
  loci_names <- paste0("locus", 1:500)
  chromosomes <- c(rep(1, 255),rep(2, 245))
  #positions <- rep(0, 500)
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
  test <- gen_tibble(x = example_genotypes, loci = example_loci, indiv_meta = example_indiv_meta, quiet = TRUE)

  # gen_tibble with no missingness runs
  test_pca <- test %>% gt_pca_autoSVD(verbose = FALSE)

  # Now try with imputed data
  library(bigsnpr)
  bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
  missing_gt <- gen_tibble(bed_file,  backingfile = tempfile("missing_"),quiet=TRUE)
  expect_error(missing_gt %>% gt_pca_autoSVD(),
                "You can't have missing")
  missing_gt <- gt_impute_simple(missing_gt, method = "mode")

  # we encounter roll_size error
  expect_error(missing_pca <- missing_gt %>% gt_pca_autoSVD(verbose = FALSE), "Parameter 'size' is too large.")

  # adjusting roll_size fixes the error
  test_pca_roll0 <- missing_gt %>% gt_pca_autoSVD(roll_size = 0, verbose = FALSE)
  expect_s3_class(test_pca_roll0, "gt_pca")
  test_pca_roll7 <- missing_gt %>% gt_pca_autoSVD(roll_size = 7, verbose = FALSE)
  expect_s3_class(test_pca_roll7, "gt_pca")


  #testing with the families dataset too
  bed_file <- system.file("extdata/related", "families.bed", package = "tidypopgen")
  families <- gen_tibble(bed_file,  backingfile = tempfile("families"),quiet=TRUE, valid_alleles = c("2","1"))
  expect_error( families %>% gt_pca_autoSVD(),
                "You can't have missing")
  families <- gt_impute_simple(families, method = "mode")

  # the same error occurs
  expect_error(missing_pca <- families %>% gt_pca_autoSVD(verbose = FALSE), "Parameter 'size' is too large.")

  # adjusting roll_size fixes the error
  test_pca_families_roll7 <- families %>% gt_pca_autoSVD(roll_size = 7, verbose = FALSE)
  expect_s3_class(test_pca_families_roll7, "gt_pca")
})



test_that("fit_gt_pca_and_predict_splitted_data",{
  bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
  missing_gt <- gen_tibble(bed_file,  backingfile = tempfile("missing_"),quiet=TRUE)
  # create a fake ancient set by subsetting
  ancient_gt <- missing_gt[1:20,]
  # now extract the modern data (to be imputed)
  modern_gt <- missing_gt[-c(1:20),]

  modern_gt <- gt_impute_simple(modern_gt, method = "mode")
  modern_pca <- modern_gt %>% gt_pca_partialSVD()
  # if we just try to predict, we find that the new data have missing data
  new_pred <- predict(modern_pca, new_data = ancient_gt, project_method = "simple")
  expect_true(all(dim(new_pred)==c(20,10)))
  # now raise an error if we don't impute to the mean
  expect_error(predict(modern_pca, new_data = ancient_gt, project_method = "none"),
               "You can't have missing values in 'X'")
  # least squares prediction
  lsq_pred <- predict(modern_pca, new_data = ancient_gt, project_method = "least_squares")
  expect_true(all(dim(lsq_pred)==c(20,2)))
})

test_that("PCA functions work with loci out of order",{
  bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
  missing_gt <- gen_tibble(bed_file,  backingfile = tempfile("missing_"),quiet=TRUE)
  missing_gt <- gt_impute_simple(missing_gt, method = "mode")
  missing_part_pca1 <- missing_gt %>% gt_pca_partialSVD()
  missing_rand_pca1 <- missing_gt %>% gt_pca_randomSVD()

  # now shuffle the loci
  loci <- missing_gt %>% show_loci()
  loci <- loci[sample(nrow(loci)),]
  show_loci(missing_gt) <- loci

  # Rerun PCA
  missing_part_pca2 <- missing_gt %>% gt_pca_partialSVD()
  missing_rand_pca2 <- missing_gt %>% gt_pca_randomSVD()

  expect_equal(missing_part_pca1$loadings, missing_part_pca2$loadings)
  expect_equal(missing_rand_pca1$loadings, missing_rand_pca2$loadings)
})

test_that("PCA computes total variance correctly",{
  bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
  missing_gt <- gen_tibble(bed_file,  backingfile = tempfile("missing_"),quiet=TRUE)
  missing_gt <- gt_impute_simple(missing_gt, method = "mode")
  missing_part_pca1 <- missing_gt %>% gt_pca_partialSVD()
  expect_true ("square_frobenious" %in% names(missing_part_pca1))
  tidy_pca_out <- tidy(missing_part_pca1,matrix = "eigenvalues")
  expect_true ("cumulative" %in% names(tidy_pca_out))
  # now repeat without estimating variance
  missing_part_pca2 <- missing_gt %>% gt_pca_partialSVD(square_frobenious = FALSE)
  expect_false ("square_frobenious" %in% names(missing_part_pca2))
  tidy_pca_out <- tidy(missing_part_pca2,matrix = "eigenvalues")
  expect_false ("cumulative" %in% names(tidy_pca_out))
  # TODO we should repeat the PCA on a toy dataset to check that it is equivalent to a full PCA done by hand

})

