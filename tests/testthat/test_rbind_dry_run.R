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
  index_pair_target <- na.omit(report$target[, c("new_id", "name")])
  index_pair_ref <- na.omit(report$ref[, c("new_id", "name")])
  index_pair_target <- index_pair_target[order(index_pair_target$new_id), ]
  index_pair_ref <- index_pair_ref[order(index_pair_ref$new_id), ]

  # check the list of new_id and name are now equal in both outputs
  expect_gt(nrow(index_pair_target), 0)
  expect_identical(index_pair_target$new_id, index_pair_ref$new_id)
  expect_identical(index_pair_target$name, index_pair_ref$name)

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

# Loci with an unobserved allele ------------------------------------------
#
# A locus can carry an unobserved allele on either side of a merge, and the
# gap can sit in either allele column: gen_tibble_bed() maps allele_ref to
# the bim file's allele2 column and allele_alt to allele1, so a missing code
# in either column of the source file can arrive in either column here.
# When the other dataset has both alleles, the gap is inferred from them,
# whichever column it is in. When both datasets have a gap, there is nothing
# to infer it from, and the locus is kept (with its gap) only if the observed
# alleles match. With flip_strand = TRUE, a kept locus whose allele_ref is
# missing in both datasets is swapped, so that the observed allele becomes
# allele_ref; with flip_strand = FALSE, the gap is left in allele_ref. A gap
# must never be carried into the report as NA (which propagates into the
# logical subscripting used downstream), nor filled with the internal missing
# code "0".
#
# Note when extending this fixture: strand-ambiguous pairs (A/T, C/G) are
# removed when flip_strand = TRUE, and that applies to the alleles *after*
# any missing allele has been resolved from the other dataset. Every pair
# below is chosen to be non-ambiguous post-resolution.

unobs_indiv_ref <- data.frame(
  id = c("a", "b", "c"), population = c("pop1", "pop1", "pop2")
)
unobs_indiv_tgt <- data.frame(
  id = c("x", "y", "z"), population = c("pop3", "pop3", "pop4")
)
unobs_geno <- rbind(
  c(1, 1, 0, 1, 1, 0, 0, 1),
  c(2, 1, 0, 0, 0, 0, 1, 2),
  c(2, 2, 0, 0, 1, 1, 2, 0)
)
unobs_base_loci <- data.frame(
  name = paste0("rs", 1:8),
  chromosome = paste0("chr", c(1, 1, 1, 1, 2, 2, 2, 2)),
  position = as.integer(c(3, 5, 65, 343, 23, 456, 500, 600)),
  genetic_dist = as.double(rep(0, 8))
)

unobs_loci_ref <- cbind(unobs_base_loci, data.frame(
  allele_ref = c("A", "T", "C", NA, "A", "T", NA, NA),
  allele_alt = c("G", "C", NA, "G", "G", "C", "G", "C")
))
unobs_loci_tgt <- cbind(unobs_base_loci, data.frame(
  allele_ref = c("A", "C", "C", "A", "G", "A", "G", NA),
  allele_alt = c("G", "T", "T", "G", NA, "G", NA, "C")
))

unobs_gt_ref <- gen_tibble(
  x = unobs_geno, loci = unobs_loci_ref, indiv_meta = unobs_indiv_ref,
  valid_alleles = c("A", "T", "C", "G"), quiet = TRUE
)
unobs_gt_tgt <- gen_tibble(
  x = unobs_geno, loci = unobs_loci_tgt, indiv_meta = unobs_indiv_tgt,
  valid_alleles = c("A", "T", "C", "G"), quiet = TRUE
)

unobs_kept_loci <- paste0("rs", 1:8)

test_that("an unobserved allele in either column does not error", {
  # the gap sits in allele_ref for rs4, rs7 and rs8 in the reference, and in
  # allele_alt for rs3 (reference), rs5 and rs7 (target)
  for (flip_strand in c(TRUE, FALSE)) {
    expect_no_error(
      unobs_report <- rbind_dry_run(unobs_gt_ref, unobs_gt_tgt,
        flip_strand = flip_strand, quiet = TRUE
      )
    )
    expect_no_error(
      rbind(unobs_gt_ref, unobs_gt_tgt,
        flip_strand = flip_strand, quiet = TRUE,
        backingfile = tempfile("test_unobserved_allele_")
      )
    )
    # unobserved alleles are never filled with the missing code, and rs7 has
    # gaps in opposite columns on the two sides, which must not be used to
    # "resolve" each other
    for (side in c("target", "ref")) {
      expect_false(any(unobs_report[[side]]$missing_allele %in% "0"))
      rs7 <- unobs_report[[side]][unobs_report[[side]]$name == "rs7", ]
      expect_true(is.na(rs7$missing_allele))
      expect_true(is.na(rs7$missing_allele_col))
    }
  }
})

test_that("an unobserved allele_alt is resolved from the other dataset", {
  unobs_report <- rbind_dry_run(unobs_gt_ref, unobs_gt_tgt,
    flip_strand = TRUE, quiet = TRUE
  )
  # rs3's gap is on the reference side, rs5's on the target side
  ref <- unobs_report$ref
  tgt <- unobs_report$target
  expect_equal(ref$missing_allele[ref$name == "rs3"], "T")
  expect_equal(ref$missing_allele_col[ref$name == "rs3"], "allele_alt")
  expect_equal(tgt$missing_allele[tgt$name == "rs5"], "A")
  expect_equal(tgt$missing_allele_col[tgt$name == "rs5"], "allele_alt")
  expect_false(is.na(tgt$new_id[tgt$name == "rs3"]))
  expect_false(is.na(tgt$new_id[tgt$name == "rs5"]))
})

test_that("an unobserved allele_ref is resolved from the other dataset", {
  # rs4 is NA/G in the reference and A/G in the target, so the missing
  # allele_ref can only be A
  for (flip_strand in c(TRUE, FALSE)) {
    unobs_report <- rbind_dry_run(unobs_gt_ref, unobs_gt_tgt,
      flip_strand = flip_strand, quiet = TRUE
    )
    ref <- unobs_report$ref
    tgt <- unobs_report$target
    expect_equal(ref$missing_allele[ref$name == "rs4"], "A")
    expect_equal(ref$missing_allele_col[ref$name == "rs4"], "allele_ref")
    expect_false(is.na(tgt$new_id[tgt$name == "rs4"]))
    expect_false(tgt$to_flip[tgt$name == "rs4"])
    expect_false(tgt$to_swap[tgt$name == "rs4"])
  }
})

test_that("swaps and flips are still detected alongside unobserved alleles", {
  # with flip_strand, rs7 and rs8 (allele_ref missing in both datasets) are
  # swapped in the reference, and the target follows its new orientation: rs7
  # (G/NA in the target) now matches it as is, and rs8 (NA/C) needs a swap
  unobs_report <- rbind_dry_run(unobs_gt_ref, unobs_gt_tgt,
    flip_strand = TRUE, quiet = TRUE
  )
  ref <- unobs_report$ref
  tgt <- unobs_report$target
  expect_setequal(ref$name[ref$to_swap], c("rs7", "rs8"))
  expect_setequal(tgt$name[tgt$to_swap], c("rs2", "rs5", "rs8"))
  expect_setequal(tgt$name[tgt$to_flip], "rs6")
  expect_setequal(tgt$name[!is.na(tgt$new_id)], unobs_kept_loci)

  # without flip_strand, nothing is swapped in the reference, and rs6 (on the
  # opposite strand) is dropped
  unobs_report <- rbind_dry_run(unobs_gt_ref, unobs_gt_tgt,
    flip_strand = FALSE, quiet = TRUE
  )
  ref <- unobs_report$ref
  tgt <- unobs_report$target
  expect_false(any(ref$to_swap))
  expect_setequal(tgt$name[tgt$to_swap], c("rs2", "rs5", "rs7"))
  expect_false(any(tgt$to_flip))
  expect_setequal(
    tgt$name[!is.na(tgt$new_id)],
    setdiff(unobs_kept_loci, "rs6")
  )
})

test_that("the summary reports the loci swapped in the reference", {
  unobs_report <- rbind_dry_run(unobs_gt_ref, unobs_gt_tgt,
    flip_strand = TRUE, quiet = TRUE
  )
  expect_output(
    summary(unobs_report),
    "2 had allele_ref missing in both datasets, and were swapped"
  )
  unobs_report <- rbind_dry_run(unobs_gt_ref, unobs_gt_tgt,
    flip_strand = FALSE, quiet = TRUE
  )
  expect_false(any(grepl(
    "allele_ref missing in both datasets",
    capture.output(summary(unobs_report))
  )))
})

test_that("rbind of loci with unobserved alleles gives a complete object", {
  # resolved gaps are filled in the right column, and the missing code "0"
  # never appears. With flip_strand, the gaps that could not be resolved (rs7
  # and rs8) are swapped into allele_alt; without it, they stay in allele_ref,
  # and rs6 (on the opposite strand) is dropped
  expected <- list(
    "TRUE" = list(
      name = unobs_kept_loci,
      allele_ref = c("A", "T", "C", "A", "A", "T", "G", "C"),
      allele_alt = c("G", "C", "T", "G", "G", "C", NA, NA)
    ),
    "FALSE" = list(
      name = setdiff(unobs_kept_loci, "rs6"),
      allele_ref = c("A", "T", "C", "A", "A", NA, NA),
      allele_alt = c("G", "C", "T", "G", "G", "G", "C")
    )
  )
  # recode the dosage (2 - g) of the loci whose alleles are swapped
  swap_geno <- function(geno, to_swap) {
    geno[, to_swap] <- 2 - geno[, to_swap]
    geno
  }
  for (flip_strand in c(TRUE, FALSE)) {
    unobs_report <- rbind_dry_run(unobs_gt_ref, unobs_gt_tgt,
      flip_strand = flip_strand, quiet = TRUE
    )
    merged <- rbind(unobs_gt_ref, unobs_gt_tgt,
      flip_strand = flip_strand, quiet = TRUE,
      backingfile = tempfile("test_unobserved_allele_")
    )
    expect_equal(nrow(merged), nrow(unobs_gt_ref) + nrow(unobs_gt_tgt))
    merged_loci <- show_loci(merged)
    expected_loci <- expected[[as.character(flip_strand)]]
    expect_equal(merged_loci$name, expected_loci$name)
    expect_equal(merged_loci$allele_ref, expected_loci$allele_ref)
    expect_equal(merged_loci$allele_alt, expected_loci$allele_alt)
    # the genotypes are recoded exactly for the loci flagged as to_swap in the
    # report (the two datasets have the same loci, in the same order)
    kept <- unobs_kept_loci %in% expected_loci$name
    expected_geno <- rbind(
      swap_geno(unobs_geno, unobs_report$ref$to_swap)[, kept],
      swap_geno(unobs_geno, unobs_report$target$to_swap)[, kept]
    )
    merged_geno <- show_genotypes(merged)
    dimnames(merged_geno) <- NULL
    expect_equal(merged_geno, expected_geno)
  }
})

test_that("a flipped locus missing allele_ref in both datasets is swapped", {
  # s1 is NA/G in the reference and NA/C in the target: the target is on the
  # opposite strand, and the observed allele becomes allele_ref in both
  # datasets. s2 (A/G vs T/C) is a plain flip, as a control
  make_gt <- function(allele_ref, allele_alt, ids) {
    gen_tibble(
      x = rbind(c(0, 1), c(1, 2), c(2, 0)),
      loci = data.frame(
        name = c("s1", "s2"),
        chromosome = "chr1",
        position = as.integer(c(1, 2)),
        genetic_dist = as.double(c(0, 0)),
        allele_ref = allele_ref,
        allele_alt = allele_alt
      ),
      indiv_meta = data.frame(id = ids, population = "pop1"),
      valid_alleles = c("A", "T", "C", "G"),
      quiet = TRUE
    )
  }
  gt_ref <- make_gt(c(NA, "A"), c("G", "G"), c("a", "b", "c"))
  gt_tgt <- make_gt(c(NA, "T"), c("C", "C"), c("x", "y", "z"))

  report <- rbind_dry_run(gt_ref, gt_tgt, flip_strand = TRUE, quiet = TRUE)
  expect_equal(report$ref$to_swap, c(TRUE, FALSE))
  expect_equal(report$target$to_flip, c(TRUE, TRUE))
  expect_equal(report$target$to_swap, c(TRUE, FALSE))

  merged <- rbind(gt_ref, gt_tgt,
    flip_strand = TRUE, quiet = TRUE,
    backingfile = tempfile("test_unobserved_allele_")
  )
  expect_equal(show_loci(merged)$allele_ref, c("G", "A"))
  expect_equal(show_loci(merged)$allele_alt, c(NA, "G"))
  merged_geno <- show_genotypes(merged)
  dimnames(merged_geno) <- NULL
  # s1 is recoded (2 - g) in both datasets, s2 in neither
  expected_geno <- cbind(c(2, 1, 0), c(1, 2, 0))
  expect_equal(merged_geno, rbind(expected_geno, expected_geno))
})


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
