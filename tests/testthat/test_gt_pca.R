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
                        predict(missing_pca, new_data = missing_gt,
                                impute_to_center = FALSE),
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
                                impute_to_center = FALSE),
                        check.attributes=FALSE))

})

# TODO we should test gt_pca_autoSVD(), as the loci have to be subset within
# the object

test_that("gt_pca_autoSVD problem ",{

  # Generate example_indiv_meta data frame
  individual_ids <- paste0("indiv", 1:200)
  populations <- sample(c("pop1", "pop2", "pop3", "pop4"), size = 200, replace = TRUE)
  example_indiv_meta <- data.frame(id = individual_ids, population = populations)

  # Generate example_genotypes matrix
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


  test <- gen_tibble(x = example_genotypes, loci = example_loci, indiv_meta = example_indiv_meta, quiet = TRUE)

  #works fine
  #test_pca <- test %>% gt_pca_autoSVD()




  # Now try with imputed data
  bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
  missing_gt <- gen_tibble(bed_file,  backingfile = tempfile("missing_"),quiet=TRUE)
  expect_error(missing_gt %>% gt_pca_partialSVD(),
                "You can't have missing")
  missing_gt <- gt_impute_simple(missing_gt, method = "mode")
  expect_error(missing_pca <- missing_gt %>% gt_pca_autoSVD(), "Parameter 'size' is too large.")

  #Issue for bigutilsr below:
  #https://github.com/privefl/bigutilsr/issues/2

  #Error in bigutilsr::rollmean(S.col[ind], roll.size) :
  #Parameter 'size' is too large.

  #the error persists when manually adjusting roll_size

  expect_error(missing_pca <- missing_gt %>% gt_pca_autoSVD(roll_size = 10), "Parameter 'size' is too large.")
  expect_error(missing_pca <- missing_gt %>% gt_pca_autoSVD(roll_size = 1000), "Parameter 'size' is too large.")

  #until adjusting roll_size to 0
  #missing_gt %>% gt_pca_autoSVD(roll_size = 0)
  #or strangely up to 7?
  #missing_gt %>% gt_pca_autoSVD(roll_size = 7)



  #loci info is the same
  #summary(test %>% show_loci())
  #summary(missing_gt %>% show_loci())

  #testing with the families dataset too
  bed_file <- system.file("extdata/related", "families.bed", package = "tidypopgen")
  missing_gt <- gen_tibble(bed_file,  backingfile = tempfile("families"),quiet=TRUE, valid_alleles = c("2","1"))
  expect_error( missing_gt %>% gt_pca_partialSVD(),
                "You can't have missing")
  missing_gt <- gt_impute_simple(missing_gt, method = "mode")

  # the same error occurs
  expect_error(missing_pca <- missing_gt %>% gt_pca_autoSVD(), "Parameter 'size' is too large.")

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
  ancient_pred <- predict(modern_pca, new_data = ancient_gt)
  expect_true(all(dim(ancient_pred)==c(20,10)))
  # now raise an error if we don't impute to the mean
  expect_error(predict(modern_pca, new_data = ancient_gt, impute_to_center = FALSE),
               "You can't have missing values in 'X'")
})
