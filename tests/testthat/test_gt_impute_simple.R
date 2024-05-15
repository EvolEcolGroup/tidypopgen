test_that("impute and use the imputation",{
  bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
  missing_gt <- gen_tibble(bed_file,  backingfile = tempfile("missing_"), quiet = TRUE)
  # we get errors because of missing values
  expect_error(missing_gt %>% gt_pca_partialSVD(),
               "You can't have missing values")
  expect_false(gt_has_imputed(missing_gt))
  expect_error(gt_uses_imputed(missing_gt),
               "this dataset does not have any imputed")
  # now impute
  missing_gt <- gt_impute_simple(missing_gt)
  # we have imputed
  expect_true(gt_has_imputed(missing_gt))
  # but don't use it by default
  expect_false(gt_uses_imputed(missing_gt))
  # now we return a pca successfully
  expect_true(inherits(missing_gt %>% gt_pca_partialSVD(),"gt_pca"))
  # simple error message
  expect_error(gt_set_imputed(missing_gt),
               "set should be either")
})


test_that("gt_impute imputes properly",{
  test_indiv_meta <- data.frame (id=c("a","b","c","d","e","f"),
                                 population = c("pop1","pop1","pop2","pop1","pop1","pop2"))

  test_genotypes <- matrix(c(
    0, 2, 1, 1, 0,  #
    0, 0, 2, 0, 1,  #
    2, 0, 0, 1, 1,  #
    1, 1, 2, 2, 2,  #
    0, 0, 2, 1, 1,  #
    NA,NA,NA,NA,NA
  ), nrow = 6, byrow = TRUE)



  test_loci <- data.frame(name=paste0("rs",1:5),
                          chromosome=c(1,1,1,2,2),
                          position=c(3,65,343,23,456),
                          genetic_dist = as.integer(rep(0,5)),
                          allele_ref = c("A","T","C","G","C"),
                          allele_alt = c("T","C", NA,"C","G"))

  test_gt <- gen_tibble(x = test_genotypes, loci = test_loci, indiv_meta = test_indiv_meta, quiet = TRUE)

  imputed_gt_mode <- gt_impute_simple(test_gt, method = "mode")
  imputed_gt_mean0 <- gt_impute_simple(test_gt, method = "mean0")
  imputed_gt_mean2 <- gt_impute_simple(test_gt, method = "mean2")
  imputed_gt_random <- gt_impute_simple(test_gt, method = "random")

  #check imputation
  expect_false(gt_has_imputed(test_gt))
  expect_true(gt_has_imputed(imputed_gt_mode))

  #test errors on non-imputed set
  expect_error(gt_uses_imputed(test_gt),"this dataset does not have any imputed values")
  expect_error(gt_set_imputed(test_gt, TRUE),"this dataset does not have imputed values")

  #set imputation
  gt_set_imputed(imputed_gt_mode, TRUE)
  gt_set_imputed(imputed_gt_mean0, TRUE)
  gt_set_imputed(imputed_gt_mean2, TRUE)
  gt_set_imputed(imputed_gt_random, TRUE)

  #test error trying to impute an already imputed set
  expect_error(gt_impute_simple(imputed_gt_mode),"object x is already imputed")

  #check there are no missing values after imputation
  expect_false(any(is.na(show_genotypes(imputed_gt_mode))))
  expect_false(any(is.na(show_genotypes(imputed_gt_mean0))))
  expect_false(any(is.na(show_genotypes(imputed_gt_mean2))))
  expect_false(any(is.na(show_genotypes(imputed_gt_random))))

  #check genotypes
  show_genotypes(imputed_gt_mode)
  show_genotypes(imputed_gt_mean0)
  show_genotypes(imputed_gt_mean2)

  #check imputed 'mean0' method
  means <- round(colMeans(test_genotypes, na.rm = TRUE), digit = 0)
  expect_true(all(means == show_genotypes(imputed_gt_mean0)[6,]))

  #Check imputed 'mode' method
  mode_function <- function(x){
    unique_x <- unique(x)
    return(unique_x[which.max(tabulate(match(x, unique_x)))])

  }

  modes <- apply(test_genotypes, 2, mode_function)
  expect_true(all(show_genotypes(imputed_gt_mode)[6,] == modes))

})



