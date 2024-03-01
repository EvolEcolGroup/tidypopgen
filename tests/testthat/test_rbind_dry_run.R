#reference file
raw_path_pop_b <- system.file("extdata/pop_b.raw", package = "tidypopgen")
map_path_pop_b <- system.file("extdata/pop_b.map", package = "tidypopgen")
pop_b_gen <- read_plink_raw(file = raw_path_pop_b, map_file = map_path_pop_b, quiet = TRUE)

#target file
raw_path_pop_a <- system.file("extdata/pop_a.raw", package = "tidypopgen")
map_path_pop_a <- system.file("extdata/pop_a.map", package = "tidypopgen")
pop_a_gen <- read_plink_raw(file = raw_path_pop_a, map_file = map_path_pop_a, quiet = TRUE)

#create merge report
report <- rbind_dry_run(pop_b_gen, pop_a_gen, flip_strand = TRUE,
                        remove_ambiguous = TRUE, quiet = TRUE)


testthat::test_that("merge report detects matching rsID's correctly",{

  #check new_id index
  #exclude NA's - those missing in either target or ref
  index_pair_target <- na.omit(report$target_gen[,c(2,3)])
  index_pair_ref <- na.omit(report$ref_gen[,c(2,3)])

  #check the list of new_id and name are now equal in both outputs
  testthat::expect_true(all(index_pair_target[,c(1,2)] == index_pair_ref[,c(1,2)]))


  # now create report directly from the files and check that it is the same as from the genlight objects
#  report_char <- em_merge_report(target = raw_path_pop_a, ref = raw_path_pop_b, quiet = TRUE)
#  testthat::expect_identical(report, report_char)

})

testthat::test_that("merge report evaluates non-matching target loci correctly",{
  #check that NA non-matching always return FALSE to_flip and to_swap

  #create list of SNPs in target_gen that are not in ref_gen
  missing_in_ref <- subset(report$target, is.na(report$target$new_id))

  #check these return false to_flip and to_swap
  testthat::expect_false(all(missing_in_ref$to_flip))
  testthat::expect_false(all(missing_in_ref$to_swap))

})


testthat::test_that("merge report detects flip alleles correctly",{
 #Based on expectations from manual inspection of the data

  #Matching strand, matching order:
  condition1 <- subset(report$target, report$target$name %in% c("rs3094315", "rs3131972", "rs1110052"))
  testthat::expect_false(all(condition1$to_flip))
  testthat::expect_false(all(condition1$to_swap))

  #Matching strand, opposite order:
  condition2 <- subset(report$target, report$target$name %in% c("rs11240777"))
  testthat::expect_false(all(condition2$to_flip))
  testthat::expect_true(all(condition2$to_swap))

})


testthat::test_that("merge report detects opposite strand alleles correctly",{

  #Opposite strand, matching order:
  condition3 <- subset(report$target, report$target$name %in% c("rs2862633","rs28569024"))
  testthat::expect_true(all(condition3$to_flip))
  testthat::expect_false(all(condition3$to_swap))

  #Opposite strand, opposite order:
  condition4 <- subset(report$target, report$target$name %in% c("rs10106770"))
  testthat::expect_true(all(condition4$to_flip))
  testthat::expect_true(all(condition4$to_swap))

})

testthat::test_that("missing cases are given the correct allele",{

  #Manually check case rs1110052: a T/G snp, coded m/t in pop_b_gen
  missing_pop_b <- subset(report$ref, name == "rs1110052")
  testthat::expect_equal(missing_pop_b$missing_allele, "t")

  #Manually check cases rs3094315 and rs11942835:
  missing_pop_a <- subset(report$target, name %in% c("rs3094315","rs11942835"))
  testthat::expect_equal(missing_pop_a$missing_allele, c("g","c"))

  #Check no missing allele given if rsID isn't in reference set

  #rs4477212 missing allele in target data, but snp not in ref data
  missing_pop_a_non_overlapping <- subset(report$target, name == "rs4477212")
  testthat::expect_true(is.na(missing_pop_a_non_overlapping$missing_allele))

})

