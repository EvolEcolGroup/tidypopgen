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
#  report_char <- rbind_dry_run(ref = raw_path_pop_b, target = raw_path_pop_a, flip_strand = TRUE,
#                               remove_ambiguous = TRUE, quiet = TRUE)
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

testthat::test_that("merge report detects ambiguous alleles correctly",{
  #Based on expectations from manual inspection of the data

  #Ambiguous found in both sets
  ambiguous_both_sets <- subset(report$target,report$target$name == "rs1240719")
  testthat::expect_true(ambiguous_both_sets$ambiguous)
  testthat::expect_true(is.na(ambiguous_both_sets$new_id))

  #Ambiguous in the target set only
  ambiguous_target_set <- subset(report$target,report$target$name == "rs307354")
  testthat::expect_true(ambiguous_both_sets$ambiguous)
  testthat::expect_true(is.na(ambiguous_both_sets$new_id))

  #Ambiguous in the ref set only
  ambiguous_ref_set <- subset(report$target,report$target$name == "rs2843130")
  testthat::expect_true(ambiguous_both_sets$ambiguous)
  testthat::expect_true(is.na(ambiguous_both_sets$new_id))


})


testthat::test_that("merge report detects flip alleles correctly",{
 #Based on expectations from manual inspection of the data

  #Matching strand, matching order:
  condition1 <- subset(report$target, report$target$name %in% c("rs3094315", "rs3131972", "rs1110052"))
  testthat::expect_true(all(condition1$to_flip == FALSE))
  testthat::expect_true(all(condition1$to_swap == FALSE))

  #Matching strand, opposite order:
  condition2 <- subset(report$target, report$target$name %in% c("rs11240777"))
  testthat::expect_false(all(condition2$to_flip))
  testthat::expect_true(all(condition2$to_swap))

})


testthat::test_that("merge report detects opposite strand alleles correctly",{

  #Opposite strand, matching order:
  condition3 <- subset(report$target, report$target$name %in% c("rs2862633","rs28569024"))
  testthat::expect_true(all(condition3$to_flip == TRUE))
  testthat::expect_true(all(condition3$to_swap == FALSE))

  #Opposite strand, opposite order:
  condition4 <- subset(report$target, report$target$name == "rs10106770")
  testthat::expect_true(all(condition4$to_flip))
  testthat::expect_true(all(condition4$to_swap))

})

testthat::test_that("missing cases are given the correct alleles",{

  #rs4477212: missing allele in target data, but snp not in ref data

  #Expect NA
  missing_pop_a_non_overlapping <- subset(report$target, name == "rs4477212")
  testthat::expect_true(is.na(missing_pop_a_non_overlapping$missing_allele))

  #rs12124819 and rs6657048: missing in target, same order, same strand
  #Target 0 A, reference G A
  #Target 0 C, reference T C

  #Expect false to_flip and to_swap
  miss_pop_a_ordered <- subset(report$target, report$target$name %in% c("rs12124819","rs6657048"))
  testthat::expect_true(miss_pop_a_ordered$missing_allele[1] == "g")
  testthat::expect_true(miss_pop_a_ordered$missing_allele[2] == "t")
  testthat::expect_true(all(miss_pop_a_ordered$to_swap == FALSE))
  testthat::expect_true(all(miss_pop_a_ordered$to_flip == FALSE))

  #rs2488991: missing in target, different order, same strand
  #Target 0 T, reference T G

  #Expect false to_flip and true to_swap
  miss_pop_a_swapped <- subset(report$target, report$target$name %in% c("rs2488991"))
  testthat::expect_true(miss_pop_a_swapped$missing_allele == "g")
  testthat::expect_false(miss_pop_a_swapped$to_flip)
  testthat::expect_true(miss_pop_a_swapped$to_swap)

  #rs5945676: missing in target, different strand, same order
  #Target 0 T, reference C A

  #Expect true to_flip and false to_swap
  miss_pop_a_flipped_swapped <- subset(report$target, report$target$name %in% c("rs5945676"))
  testthat::expect_true(miss_pop_a_flipped_swapped$missing_allele == "g")
  testthat::expect_false(miss_pop_a_flipped_swapped$to_swap)
  testthat::expect_true(miss_pop_a_flipped_swapped$to_flip)

})

