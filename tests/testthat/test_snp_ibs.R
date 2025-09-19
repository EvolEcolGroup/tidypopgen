# create file
test_indiv_meta <- data.frame(
  id = c("a", "b", "c"),
  population = c("pop1", "pop1", "pop2")
)
test_genotypes <- rbind(
  c(1, 1, 0, 1, 1, 0),
  c(2, 1, 0, 0, 0, 0),
  c(2, 2, 0, 0, 1, 1)
)
test_loci <- data.frame(
  name = paste0("rs", 1:6),
  chromosome = paste0("chr", c(1, 1, 1, 1, 2, 2)),
  position = as.integer(c(3, 5, 65, 343, 23, 456)),
  genetic_dist = as.double(rep(0, 6)),
  allele_ref = c("A", "T", "C", "G", "C", "T"),
  allele_alt = c("T", "C", NA, "C", "G", "A")
)

test_gt <- gen_tibble(
  x = test_genotypes,
  loci = test_loci,
  indiv_meta = test_indiv_meta,
  quiet = TRUE
)


test_that("snp_ibs and pairwise_ibs computes ibs correctly", {
  test_fbm <- tidypopgen:::.gt_get_fbm(test_gt)
  test_ibs <- snp_ibs(test_fbm, type = "raw_counts")
  # compare indiv 1 vs 2
  in_common <- sum(c(1, 2, 2, 1, 1, 2))
  expect_identical(in_common, test_ibs$ibs[1, 2])
  # check that we get the same result if we split the operation into two blocks
  test_ibs_2blocks <- snp_ibs(test_fbm, block.size = 3, type = "raw_counts")
  expect_identical(test_ibs_2blocks$ibs[], test_ibs$ibs[])

  # now estimate it with gen_tibble
  test_ibs_gt <- pairwise_ibs(test_gt, type = "raw_counts")
  expect_true(all.equal(
    test_ibs$ibs[],
    test_ibs_gt$ibs[],
    check.attributes = FALSE
  ))
  test_ibs_gt_prop <- pairwise_ibs(test_gt)
  # expect_true (all.equal(test_ibs_gt_prop, test_ibs$ibs[]/test_ibs$valid_n[],
  #                                          check.attributes = FALSE))

  # use raw_counts matrices to calculate proportion
  by_hand <- as.data.frame(test_ibs$ibs[] / test_ibs$valid_n[])
  colnames(by_hand) <- c("a", "b", "c")

  # check against proportion value
  expect_true(all.equal(test_ibs_gt_prop$value[1], by_hand$a[2])) # a and b
  expect_true(all.equal(test_ibs_gt_prop$value[2], by_hand$a[3])) # a and c
  expect_true(all.equal(test_ibs_gt_prop$value[3], by_hand$b[3])) # b and c

  # now subset to the first and second individual, and a subset of loci
  test_gt_sub <- test_gt[c(1, 3), ]
  loci_subset <- c(2, 3, 5, 6)
  in_common_1vs3 <- sum(c(1, 1, 2, 1, 2, 1)[loci_subset])

  test_gt_sub <- test_gt_sub %>% select_loci(dplyr::all_of(loci_subset))
  test_ibs_sub <- pairwise_ibs(test_gt_sub, type = "raw_counts")
  expect_identical(in_common_1vs3, test_ibs_sub$ibs[1, 2])
})


test_that("snp_ibs as.counts = FALSE gives the same results as plink", {
  # Read in results from:
  # plink --bfile families --distance square flat-missing ibs
  plink_ibs <- read.table(
    system.file("extdata/related/test_plinkIBS.mibs", package = "tidypopgen"),
    header = FALSE
  )
  # Transform to matrix
  plink_matrix <- unname(as.matrix(plink_ibs))

  # Create gentibble for the same data
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

  # Get snp_ibs results
  families_fbm <- tidypopgen:::.gt_get_fbm(families)
  tidy_ibs <- snp_ibs(families_fbm)

  # Check both are numeric and round
  tidy_ibs <- as.numeric(tidy_ibs)
  plink_matrix <- as.numeric(plink_matrix)
  tidy_ibs <- round(tidy_ibs, 6)

  # Check matrices are equal
  expect_true(all.equal(tidy_ibs, plink_matrix))
})
