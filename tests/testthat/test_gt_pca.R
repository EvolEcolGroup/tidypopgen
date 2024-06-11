test_that("fit_gt_pca_and_predict",{
  bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
  missing_gt <- gen_tibble(bed_file,  backingfile = tempfile("missing_"),quiet=TRUE)
  expect_error( missing_gt %>% gt_pca_partialSVD(),
                "You can't have missing")
  missing_gt <- gt_impute_simple(missing_gt, method = "mode")
  missing_pca <- missing_gt %>% gt_pca_partialSVD()
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
                        predict(missing_sub_pca, new_data = missing_gt),
                        check.attributes=FALSE))

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
  expect_error(predict(modern_pca, new_data = ancient_gt),
                        "You can't have")
})
