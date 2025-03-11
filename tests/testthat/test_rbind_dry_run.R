# reference file
raw_path_pop_b <- system.file("extdata/pop_b.bed", package = "tidypopgen")
bigsnp_path_b <- bigsnpr::snp_readBed(
  raw_path_pop_b,
  backingfile = tempfile("test_b_")
)
pop_b_gt <- gen_tibble(bigsnp_path_b, quiet = TRUE)
# target file
raw_path_pop_a <- system.file("extdata/pop_a.bed", package = "tidypopgen")
bigsnp_path_a <- bigsnpr::snp_readBed(
  raw_path_pop_a,
  backingfile = tempfile("test_a_")
)
pop_a_gt <- gen_tibble(bigsnp_path_a, quiet = TRUE)


# create merge report
report <- rbind_dry_run(pop_b_gt, pop_a_gt, flip_strand = TRUE, quiet = TRUE)


test_that("merge report detects matching rsID's correctly", {
  # check new_id index
  # exclude NA's - those missing in either target or ref
  index_pair_target <- na.omit(report$target_gen[, c(2, 3)])
  index_pair_ref <- na.omit(report$ref_gen[, c(2, 3)])

  # check the list of new_id and name are now equal in both outputs
  expect_true(all(
    index_pair_target[, c(1, 2)] == index_pair_ref[, c(1, 2)]
  ))

  # now create report directly from the bim files and check that it is the same
  # as from the gen_tibble objects
  #  report_char <- rbind_dry_run(ref = raw_path_pop_b, #nolint start
  #                               target = raw_path_pop_a,
  #                               flip_strand = TRUE,
  #                               quiet = TRUE)
  #  expect_identical(report, report_char) #nolint end
})

test_that("merge report evaluates non-matching target loci correctly", {
  # check that NA non-matching always return FALSE to_flip and to_swap

  # create list of SNPs in target_gen that are not in ref_gen
  missing_in_ref <- subset(report$target, is.na(report$target$new_id))

  # check these return false to_flip and to_swap
  expect_false(all(missing_in_ref$to_flip))
  expect_false(all(missing_in_ref$to_swap))
})

test_that("merge report detects ambiguous alleles correctly", {
  # Based on expectations from manual inspection of the data

  # Ambiguous found in both sets
  ambiguous_both_sets <- subset(
    report$target,
    report$target$name == "rs1240719"
  )
  expect_true(ambiguous_both_sets$ambiguous)
  expect_true(is.na(ambiguous_both_sets$new_id))

  # Ambiguous in the target set only
  ambiguous_target_set <- subset(
    report$target,
    report$target$name == "rs307354"
  )
  expect_true(ambiguous_both_sets$ambiguous)
  expect_true(is.na(ambiguous_both_sets$new_id))

  # Ambiguous in the ref set only
  ambiguous_ref_set <- subset(report$target, report$target$name == "rs2843130")
  expect_true(ambiguous_both_sets$ambiguous)
  expect_true(is.na(ambiguous_both_sets$new_id))
})


test_that("merge report detects flip alleles correctly", {
  # Based on expectations from manual inspection of the data

  # Matching strand, matching order:
  condition1 <- subset(
    report$target,
    report$target$name %in% c("rs3094315", "rs3131972", "rs1110052")
  )
  expect_true(all(condition1$to_flip == FALSE))
  expect_true(all(condition1$to_swap == FALSE))

  # Matching strand, opposite order:
  condition2 <- subset(report$target, report$target$name %in% c("rs11240777"))
  expect_false(all(condition2$to_flip))
  expect_true(all(condition2$to_swap))
})


test_that("merge report detects opposite strand alleles correctly", {
  # Opposite strand, matching order:
  condition3 <- subset(
    report$target,
    report$target$name %in% c("rs2862633", "rs28569024")
  )
  expect_true(all(condition3$to_flip == TRUE))
  expect_true(all(condition3$to_swap == FALSE))

  # Opposite strand, opposite order:
  condition4 <- subset(report$target, report$target$name == "rs10106770")
  expect_true(all(condition4$to_flip))
  expect_true(all(condition4$to_swap))
})

test_that("missing cases are given the correct alleles", {
  # rs4477212: missing allele in target data, but snp not in ref data

  # Expect NA
  missing_pop_a_non_overlapping <- subset(report$target, name == "rs4477212")
  expect_true(is.na(missing_pop_a_non_overlapping$missing_allele))

  # rs12124819 and rs6657048: missing in target, same order, same strand
  # Target 0 A, reference G A
  # Target 0 C, reference T C

  # Expect false to_flip and to_swap
  miss_pop_a_ordered <- subset(
    report$target,
    report$target$name %in% c("rs12124819", "rs6657048")
  )
  expect_true(miss_pop_a_ordered$missing_allele[1] == "G")
  expect_true(miss_pop_a_ordered$missing_allele[2] == "T")
  expect_true(all(miss_pop_a_ordered$to_swap == FALSE))
  expect_true(all(miss_pop_a_ordered$to_flip == FALSE))

  # rs2488991: missing in target, different order, same strand
  # Target 0 T, reference T G

  # Expect false to_flip and true to_swap
  miss_pop_a_swapped <- subset(
    report$target,
    report$target$name %in% c("rs2488991")
  )
  expect_true(miss_pop_a_swapped$missing_allele == "G")
  expect_false(miss_pop_a_swapped$to_flip)
  expect_true(miss_pop_a_swapped$to_swap)

  # rs5945676: missing in target, different strand, same order
  # Target 0 T, reference C A

  # Expect true to_flip and false to_swap
  miss_pop_a_flipped_swapped <- subset(
    report$target,
    report$target$name %in% c("rs5945676")
  )
  expect_true(miss_pop_a_flipped_swapped$missing_allele == "G")
  expect_false(miss_pop_a_flipped_swapped$to_swap)
  expect_true(miss_pop_a_flipped_swapped$to_flip)
})

#
#
# #reference file reordered
# raw_path_reordered_pop_b <- #nolint start
#       system.file("extdata/pop_b_reordered.raw", package = "tidypopgen")
# map_path_reordered_pop_b <-
#       system.file("extdata/pop_b_reordered.map", package = "tidypopgen")
# pop_b_gen_reordered <-
#       read_plink_raw(file = raw_path_reordered_pop_b,
#                      map_file = map_path_reordered_pop_b, quiet = TRUE) #nolint end
#
#
# test_that("reordering",{ #nolint start
#
#   #create merge report
#   report <- rbind_dry_run(pop_b_gt, pop_a_gt, flip_strand = TRUE,
#                           quiet = TRUE)
#
#   #create a new merge report with a dataset in a different order
#   report2 <- rbind_dry_run(pop_b_gen_reordered, pop_a_gt, flip_strand = TRUE,
#                           quiet = TRUE)
#
#   #Store the results of the merge report for the target data
#   report_original <- report$target
#   report_new_order <- report2$target
#
#   #Order both reports
#   report_original <- report_original[order(report_original$name),]
#   report_new_order <- report_new_order[order(report_new_order$name),]
#
#   #Deselect the new_id column
#   report_original <- report_original[c(1,3:7)]
#   report_new_order <- report_new_order[c(1,3:7)]
#
#   #Check whether merge report is the same
#   expect_identical(report_original,report_new_order)
#   #This is where the test fails
#
#
#
# }) #nolint end
