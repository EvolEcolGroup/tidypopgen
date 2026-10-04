# Loci with an unobserved allele.
#
# rbind_dry_run() replaces NA alleles with "0" before matching, because NA does
# not survive the logical subscripting used downstream. That substitution
# originally covered allele_alt only, so a missing allele arriving in allele_ref
# -- which is where gen_tibble_bed() puts the bim file's allele2 column --
# left an NA in place and the call failed under flip_strand = TRUE with
# "NAs are not allowed in subscripted assignments".
#
# Note when extending this fixture: strand-ambiguous pairs (A/T, C/G) are
# removed when flip_strand = TRUE, and that applies to the alleles *after* any
# missing allele has been resolved from the other dataset. Every pair below is
# chosen to be non-ambiguous post-resolution.

indiv_ref <- data.frame(
  id = c("a", "b", "c"), population = c("pop1", "pop1", "pop2")
)
indiv_tgt <- data.frame(
  id = c("x", "y", "z"), population = c("pop3", "pop3", "pop4")
)
geno <- rbind(
  c(1, 1, 0, 1, 1, 0),
  c(2, 1, 0, 0, 0, 0),
  c(2, 2, 0, 0, 1, 1)
)
base_loci <- data.frame(
  name = paste0("rs", 1:6),
  chromosome = paste0("chr", c(1, 1, 1, 1, 2, 2)),
  position = as.integer(c(3, 5, 65, 343, 23, 456)),
  genetic_dist = as.double(rep(0, 6))
)

#        ref        target     expected
# rs1    A/G        A/G        kept, matches as is
# rs2    T/C        C/T        kept, needs swap
# rs3    C/NA       C/T        kept, ref allele_alt resolved from target
# rs4    NA/G       A/G        dropped: gap is in allele_ref, not resolvable
# rs5    A/G        G/NA       kept, target allele_alt resolved, needs swap
# rs6    T/C        A/G        kept, needs strand flip
loci_ref <- cbind(base_loci, data.frame(
  allele_ref = c("A", "T", "C", NA, "A", "T"),
  allele_alt = c("G", "C", NA, "G", "G", "C")
))
loci_tgt <- cbind(base_loci, data.frame(
  allele_ref = c("A", "C", "C", "A", "G", "A"),
  allele_alt = c("G", "T", "T", "G", NA, "G")
))

gt_ref <- gen_tibble(
  x = geno, loci = loci_ref, indiv_meta = indiv_ref,
  valid_alleles = c("A", "T", "C", "G"), quiet = TRUE
)
gt_tgt <- gen_tibble(
  x = geno, loci = loci_tgt, indiv_meta = indiv_tgt,
  valid_alleles = c("A", "T", "C", "G"), quiet = TRUE
)

kept_loci <- c("rs1", "rs2", "rs3", "rs5", "rs6")


test_that("rbind_dry_run does not error on a missing allele in allele_ref", {
  # the regression: before the fix, the NA left in allele_ref reached a logical
  # subscript and the call failed
  expect_no_error(
    rbind_dry_run(gt_ref, gt_tgt, flip_strand = TRUE, quiet = TRUE)
  )
  expect_no_error(
    rbind_dry_run(gt_ref, gt_tgt, flip_strand = FALSE, quiet = TRUE)
  )
})

test_that("no NA leaks into the merge report", {
  report <- rbind_dry_run(gt_ref, gt_tgt, flip_strand = TRUE, quiet = TRUE)
  expect_false(anyNA(report$target$to_flip))
  expect_false(anyNA(report$target$to_swap))
  expect_false(anyNA(report$target$name))
  expect_false(anyNA(report$ref$name))
})

test_that("a missing allele_alt is still resolved from the other dataset", {
  report <- rbind_dry_run(gt_ref, gt_tgt, flip_strand = TRUE, quiet = TRUE)
  # rs3's gap is on the reference side, rs5's on the target side
  expect_equal(report$ref$missing_allele[report$ref$name == "rs3"], "T")
  expect_equal(report$target$missing_allele[report$target$name == "rs5"], "A")
  expect_false(is.na(report$target$new_id[report$target$name == "rs3"]))
  expect_false(is.na(report$target$new_id[report$target$name == "rs5"]))
})

test_that("a missing allele_ref is dropped rather than resolved", {
  # resolve_missing_alleles() only inspects allele_1 (i.e. allele_alt), so a gap
  # in allele_ref cannot be recovered. It must be dropped quietly, not error.
  report <- rbind_dry_run(gt_ref, gt_tgt, flip_strand = TRUE, quiet = TRUE)
  expect_true(is.na(report$target$new_id[report$target$name == "rs4"]))
  expect_false(report$target$to_flip[report$target$name == "rs4"])
  expect_false(report$target$to_swap[report$target$name == "rs4"])
})

test_that("normal harmonisation is unaffected", {
  report <- rbind_dry_run(gt_ref, gt_tgt, flip_strand = TRUE, quiet = TRUE)
  expect_true(report$target$to_swap[report$target$name == "rs2"])
  expect_true(report$target$to_swap[report$target$name == "rs5"])
  expect_true(report$target$to_flip[report$target$name == "rs6"])
  expect_setequal(report$target$name[!is.na(report$target$new_id)], kept_loci)
})

test_that("rbind produces a merged object with no missing alleles", {
  merged <- rbind(gt_ref, gt_tgt,
    flip_strand = TRUE, quiet = TRUE,
    backingfile = tempfile("test_missing_allele_")
  )
  expect_equal(nrow(merged), nrow(gt_ref) + nrow(gt_tgt))
  expect_setequal(loci_names(merged), kept_loci)
  expect_false(anyNA(show_loci(merged)$allele_ref))
  expect_false(anyNA(show_loci(merged)$allele_alt))
})
