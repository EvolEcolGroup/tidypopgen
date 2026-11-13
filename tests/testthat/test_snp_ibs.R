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

  # adjusted counts are the number of alleles that differ (no genotypes are
  # missing, so they are not rescaled)
  test_dist <- snp_ibs(test_fbm, type = "adjusted_counts")
  differ_by_hand <- outer(1:3, 1:3, Vectorize(function(i, j) {
    sum(abs(test_genotypes[i, ] - test_genotypes[j, ]))
  }))
  expect_equal(test_dist, differ_by_hand, ignore_attr = TRUE)
  # a and b differ by one allele at loci 1, 4 and 5
  expect_identical(test_dist[1, 2], 3)
  expect_identical(diag(test_dist), c(0, 0, 0))
  test_dist_gt <- pairwise_ibs(test_gt, type = "adjusted_counts")
  expect_identical(test_dist_gt$value, test_dist[upper.tri(test_dist)])

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

test_that("snp_ibs adjusted_counts gives the same results as plink", {
  # Read in results from (PLINK v1.9.0 64-bit, 27 Sep 2026):
  # plink --bfile families --distance square flat-missing
  # (the allele count distance, with genotypes missing in 352 of 961 loci)
  plink_dist <- read.table(
    system.file("extdata/related/test_plinkDist.dist", package = "tidypopgen"),
    header = FALSE
  )
  plink_dist <- unname(as.matrix(plink_dist))

  bed_path <- system.file(
    "extdata/related/families.bed",
    package = "tidypopgen"
  )
  families <- gen_tibble(
    bed_path,
    quiet = TRUE,
    valid_alleles = c("1", "2"),
    backingfile = tempfile()
  )
  families_fbm <- tidypopgen:::.gt_get_fbm(families)
  tidy_dist <- snp_ibs(families_fbm, type = "adjusted_counts")
  # PLINK writes 6 significant digits
  expect_equal(unname(tidy_dist), plink_dist, tolerance = 1e-5)
  expect_identical(diag(tidy_dist), rep(0, nrow(families)))
  # it is a dissimilarity in alleles: 2 * n_loci * (1 - proportion)
  expect_equal(
    tidy_dist,
    2 * count_loci(families) * (1 - snp_ibs(families_fbm))
  )
  # and pairwise_ibs gives the same values
  families_dist <- pairwise_ibs(families,
    type = "adjusted_counts",
    as_matrix = TRUE
  )
  expect_equal(unclass(families_dist), tidy_dist, ignore_attr = TRUE)
})

test_that("snp_ibs is correct across bit-packing boundaries", {
  # more than 64 individuals and loci, with missing data, so that genotypes
  # are packed into several words and blocks of individuals
  set.seed(123)
  geno_vals <- sample(c(0:2, NA), 70 * 150,
    replace = TRUE,
    prob = c(0.4, 0.3, 0.25, 0.05)
  )
  geno_mat <- matrix(geno_vals, nrow = 70)
  big_gt <- gen_tibble(
    x = geno_mat,
    loci = data.frame(
      name = paste0("rs", 1:150), chromosome = 1L, position = 1:150,
      genetic_dist = 0, allele_ref = "A", allele_alt = "T"
    ),
    indiv_meta = data.frame(id = paste0("i", 1:70), population = "p1"),
    quiet = TRUE
  )
  test_fbm <- tidypopgen:::.gt_get_fbm(big_gt)
  # reordered subsets of individuals and loci, and multiple cores
  rows <- c(70:3, 1)
  cols <- c(150:80, 2:60)
  test_ibs <- snp_ibs(test_fbm,
    ind.row = rows, ind.col = cols,
    type = "raw_counts", block.size = 70, n_cores = 2
  )
  # reference: the sum of (2 - |g_i - g_j|) over the loci valid in both
  # individuals, and 2 * the number of those loci
  g <- geno_mat[rows, cols]
  valid <- !is.na(g) + 0
  g[is.na(g)] <- 0
  abs_diff <- sapply(seq_len(nrow(g)), function(j) {
    rowSums(abs(sweep(g, 2, g[j, ])) * sweep(valid, 2, valid[j, ], "*"))
  })
  expect_identical(test_ibs$ibs[], 2 * tcrossprod(valid) - abs_diff)
  expect_identical(test_ibs$valid_n[], 2 * tcrossprod(valid))
})

test_that("pairwise_ibs uses imputed genotypes only if set to", {
  geno_na <- test_genotypes
  geno_na[2, 4] <- NA
  geno_na[3, 2] <- NA
  gt_na <- gen_tibble(
    x = geno_na, loci = test_loci, indiv_meta = test_indiv_meta,
    quiet = TRUE
  )
  ibs_na <- pairwise_ibs(gt_na, type = "raw_counts")
  gt_imputed <- gt_impute_simple(gt_na, method = "mode")
  # imputed genotypes are missing unless the gen_tibble is set to use them
  gt_set_imputed(gt_imputed, set = FALSE)
  ibs_not_used <- pairwise_ibs(gt_imputed, type = "raw_counts")
  expect_identical(ibs_not_used$ibs[], ibs_na$ibs[])
  expect_identical(ibs_not_used$valid_n[], ibs_na$valid_n[])
  # when used, they are the same as genotypes that were never missing
  gt_set_imputed(gt_imputed, set = TRUE)
  gt_filled <- gen_tibble(
    x = show_genotypes(gt_imputed), loci = test_loci,
    indiv_meta = test_indiv_meta, quiet = TRUE
  )
  ibs_used <- pairwise_ibs(gt_imputed, type = "raw_counts")
  ibs_filled <- pairwise_ibs(gt_filled, type = "raw_counts")
  expect_identical(ibs_used$ibs[], ibs_filled$ibs[])
  expect_identical(ibs_used$valid_n[], ibs_filled$valid_n[])
})

test_that("missing_correction = 'allele_freq' gives the same values as plink", {
  # Read in results from (PLINK v1.9.0 64-bit, 27 Sep 2026; see
  # data-raw/validate/create_plink_ibs_missing_data.R), for 24 individuals and
  # 300 loci with 33% of missing genotypes:
  # plink --bfile ibs_missing --distance square bin allele-ct ibs
  # [flat-missing], and, with the frequencies of the first 16 individuals
  # (the founders of ibs_missing_nonfounders.fam):
  # plink --bfile ibs_missing --fam ibs_missing_nonfounders.fam
  #   --distance square bin allele-ct
  # The .bin files hold the matrices as doubles, at full precision
  read_plink_bin <- function(file) {
    file <- system.file("extdata/related", file, package = "tidypopgen")
    matrix(readBin(file, "double", n = 24 * 24), 24, 24)
  }
  bed_path <- system.file(
    "extdata/related/ibs_missing.bed",
    package = "tidypopgen"
  )
  ibs_gt <- gen_tibble(
    bed_path,
    quiet = TRUE,
    valid_alleles = c("A", "T"),
    backingfile = tempfile()
  )
  ibs_fbm <- tidypopgen:::.gt_get_fbm(ibs_gt)
  # the weighted correction: the operations are those of PLINK, so the values
  # are identical (with several blocks and cores too)
  w_dist <- snp_ibs(ibs_fbm,
    type = "adjusted_counts",
    missing_correction = "allele_freq", block.size = 70, n_cores = 2
  )
  expect_identical(w_dist, read_plink_bin("ibs_missing.dist.bin"))
  w_prop <- pairwise_ibs(ibs_gt,
    missing_correction = "allele_freq",
    as_matrix = TRUE
  )
  expect_identical(
    unclass(unname(w_prop)),
    read_plink_bin("ibs_missing.mibs.bin")
  )
  # allele frequencies of the founders only
  founder_freqs <- loci_alt_freq(ibs_gt[1:16, ])
  w_dist_founders <- pairwise_ibs(ibs_gt,
    type = "adjusted_counts",
    missing_correction = "allele_freq", allele_freqs = founder_freqs,
    as_matrix = TRUE
  )
  expect_identical(
    unclass(unname(w_dist_founders)),
    read_plink_bin("ibs_missing_nonfounders.dist.bin")
  )
  # the flat correction (the default), up to rounding
  expect_equal(
    unname(snp_ibs(ibs_fbm, type = "adjusted_counts")),
    read_plink_bin("ibs_missing_flat.dist.bin")
  )
  expect_equal(
    unname(snp_ibs(ibs_fbm)),
    read_plink_bin("ibs_missing_flat.mibs.bin")
  )
})

test_that("missing_correction = 'allele_freq' without missing genotypes", {
  geno_na <- test_genotypes
  geno_na[2, 4] <- NA
  geno_na[3, 2] <- NA
  gt_na <- gen_tibble(
    x = geno_na, loci = test_loci, indiv_meta = test_indiv_meta,
    quiet = TRUE
  )
  # the two corrections differ with missing genotypes
  expect_false(isTRUE(all.equal(
    pairwise_ibs(gt_na, missing_correction = "allele_freq"),
    pairwise_ibs(gt_na)
  )))
  # but not without them (here, with imputed genotypes)
  gt_imputed <- gt_impute_simple(gt_na, method = "mode")
  gt_set_imputed(gt_imputed, set = TRUE)
  expect_equal(
    pairwise_ibs(gt_imputed,
      type = "adjusted_counts",
      missing_correction = "allele_freq"
    ),
    pairwise_ibs(gt_imputed, type = "adjusted_counts")
  )
  # invalid combinations of arguments
  expect_error(
    pairwise_ibs(gt_na,
      type = "raw_counts",
      missing_correction = "allele_freq"
    ),
    "can not be used with"
  )
  expect_error(
    pairwise_ibs(gt_na, allele_freqs = rep(0.5, 6)),
    "only used with"
  )
  expect_error(
    pairwise_ibs(gt_na,
      missing_correction = "allele_freq",
      allele_freqs = rep(0, 6)
    ),
    "all loci are monomorphic"
  )
})

test_that("pairwise functions delete their temporary FBMs", {
  # their temporary .bk files used to accumulate in tempdir(), where a new
  # tempfile() name could clash with them (see with_private_fbm_dir())
  test_fbm <- tidypopgen:::.gt_get_fbm(test_gt)
  n_bk <- function() length(list.files(tempdir(), pattern = "\\.bk$"))
  bk_before <- n_bk()
  snp_ibs(test_fbm, type = "proportion")
  snp_king(test_fbm)
  snp_allele_sharing(test_fbm)
  expect_identical(n_bk(), bk_before)
  # the raw counts are returned as FBMs, so their files are kept
  test_ibs <- snp_ibs(test_fbm, type = "raw_counts")
  expect_true(file.exists(test_ibs$ibs$backingfile))
  expect_true(file.exists(test_ibs$valid_n$backingfile))
})
