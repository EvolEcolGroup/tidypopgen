# the PLINK fixture: 60 individuals, two chromosomes of 100 loci in LD, with
# missing genotypes (see data-raw/validate/create_plink_ld_test_data.R)
plink_dir <- test_path("testdata", "plink_ld")
ld_gt <- gen_tibble(file.path(plink_dir, "ld.bed"),
  backingfile = tempfile("ld_"), quiet = TRUE
)

test_that("r2 over the individuals genotyped at both loci matches PLINK", {
  plink_r2 <- utils::read.table(file.path(plink_dir, "ld.ld"), header = TRUE)
  pairs <- rbind(
    match(plink_r2$SNP_A, loci_names(ld_gt)),
    match(plink_r2$SNP_B, loci_names(ld_gt))
  )
  r2 <- ld_r2_pairs_cpp(
    .gt_get_fbm(ld_gt), .gt_fbm_rows(ld_gt), .gt_fbm_cols(ld_gt), pairs
  )
  # (PLINK writes 6 digits)
  expect_equal(r2$r2, plink_r2$R2, tolerance = 1e-5)
  # the same as cor() over the individuals genotyped at both loci
  geno <- show_genotypes(ld_gt)
  shared <- !is.na(geno[, pairs[1, 1]]) & !is.na(geno[, pairs[2, 1]])
  expect_equal(
    r2$r2[1],
    stats::cor(geno[shared, pairs[1, 1]], geno[shared, pairs[2, 1]])^2
  )
  expect_identical(r2$n_shared[1], sum(shared))
})

test_that("r2 of complete and incomplete loci matches cor()", {
  # example-missing (bigsnpr): 200 individuals and 500 loci, 133 of them
  # complete; the r2 of two complete loci takes a shortcut (only the sum of
  # the products is counted), which must not be taken when only one of the
  # two loci is complete
  bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
  missing_gt <- gen_tibble(
    bed_file,
    backingfile = tempfile("missing_"),
    quiet = TRUE
  )
  geno <- show_genotypes(missing_gt)
  complete <- colSums(is.na(geno)) == 0
  # each locus with the next five
  pairs <- rbind(rep(1:495, each = 5), rep(1:495, each = 5) + 1:5)
  # pairs of two complete loci, of a complete and an incomplete locus (in
  # both orders), and of two incomplete loci are all tested
  expect_true(all(table(complete[pairs[1, ]], complete[pairs[2, ]]) > 0))
  r2 <- ld_r2_pairs_cpp(
    .gt_get_fbm(missing_gt), .gt_fbm_rows(missing_gt),
    .gt_fbm_cols(missing_gt), pairs
  )
  # the reference: cor() over the individuals genotyped at both loci
  ref <- apply(pairs, 2, function(p) {
    shared <- !is.na(geno[, p[1]]) & !is.na(geno[, p[2]])
    c(sum(shared), stats::cor(geno[shared, p[1]], geno[shared, p[2]])^2)
  })
  expect_identical(r2$n_shared, as.integer(ref[1, ]))
  expect_equal(r2$r2, ref[2, ])
})

test_that("loci_ld_prune matches plink --indep-pairwise", {
  # file name, then window_size, step_size, thr_r2 and window_in_kb
  settings <- list(
    list("w10_s1_r05", 10, 1, 0.5, FALSE),
    list("w20_s5_r02", 20, 5, 0.2, FALSE),
    list("w2_s1_r08", 2, 1, 0.8, FALSE),
    list("w50kb_s1_r02", 50, 1, 0.2, TRUE)
  )
  for (s in settings) {
    keep <- loci_ld_prune(ld_gt,
      window_size = s[[2]], step_size = s[[3]], thr_r2 = s[[4]],
      window_in_kb = s[[5]]
    )
    plink_in <- readLines(file.path(plink_dir, paste0(s[[1]], ".prune.in")))
    expect_identical(loci_names(ld_gt)[keep], plink_in)
  }
  # the indices of the loci kept
  expect_identical(
    loci_ld_prune(ld_gt,
      window_size = 10, step_size = 1, thr_r2 = 0.5, return_id = TRUE
    ),
    which(loci_ld_prune(ld_gt, window_size = 10, step_size = 1, thr_r2 = 0.5))
  )
})

test_that("r2 and pruning match PLINK with several words per locus", {
  # the individuals of the PLINK fixture stacked three times: each bit plane
  # of a locus then takes three 64-bit words (180 individuals, the last word
  # partly filled), while r2, the minor allele frequencies, and so the loci
  # kept, are those of the fixture
  geno <- show_genotypes(ld_gt)
  stacked_gt <- gen_tibble(
    x = rbind(geno, geno, geno),
    indiv_meta = data.frame(id = paste0("ind", 1:180), population = "pop"),
    loci = show_loci(ld_gt) %>% select(-big_index),
    backingfile = tempfile("ld_stacked_"),
    quiet = TRUE
  )
  # r2 of the pairs written by plink --r2
  plink_r2 <- utils::read.table(file.path(plink_dir, "ld.ld"), header = TRUE)
  pairs <- rbind(
    match(plink_r2$SNP_A, loci_names(stacked_gt)),
    match(plink_r2$SNP_B, loci_names(stacked_gt))
  )
  r2 <- ld_r2_pairs_cpp(
    .gt_get_fbm(stacked_gt), .gt_fbm_rows(stacked_gt),
    .gt_fbm_cols(stacked_gt), pairs
  )
  expect_equal(r2$r2, plink_r2$R2, tolerance = 1e-5)
  # the loci kept by plink --indep-pairwise, with windows in loci and in kb
  keep <- loci_ld_prune(stacked_gt,
    window_size = 20, step_size = 5, thr_r2 = 0.2
  )
  plink_in <- readLines(file.path(plink_dir, "w20_s5_r02.prune.in"))
  expect_identical(loci_names(stacked_gt)[keep], plink_in)
  keep <- loci_ld_prune(stacked_gt,
    window_size = 50, step_size = 1, thr_r2 = 0.2, window_in_kb = TRUE
  )
  plink_in <- readLines(file.path(plink_dir, "w50kb_s1_r02.prune.in"))
  expect_identical(loci_names(stacked_gt)[keep], plink_in)
})

test_that("loci_ld_prune checks its parameters as PLINK does", {
  expect_error(
    loci_ld_prune(ld_gt, window_size = 1, step_size = 1, thr_r2 = 0.5),
    "'window_size' should be an integer of at least 2"
  )
  expect_error(
    loci_ld_prune(ld_gt, window_size = 5, step_size = 10, thr_r2 = 0.5),
    "'step_size' should not be larger than 'window_size'"
  )
  expect_warning(
    loci_ld_prune(ld_gt,
      window_size = 50, step_size = 2, thr_r2 = 0.5, window_in_kb = TRUE
    ),
    "'step_size' should be 1 when the window is in kb"
  )
  expect_error(
    loci_ld_prune(ld_gt, window_size = 10, step_size = 1, thr_r2 = 1),
    "'thr_r2' should be at least 0, and less than 1"
  )
  expect_error(
    loci_ld_prune(ld_gt,
      window_size = 10, step_size = 1, thr_r2 = 0.5, min_shared = 1
    ),
    "'min_shared' should be an integer of at least 2"
  )
})

test_that("loci_ld_prune treats pairs with few shared individuals as asked", {
  # with a minimum above the number of individuals, no pair has an r2: with
  # "keep", only the loci with no variance are pruned (as in PLINK: the
  # monomorphic locus, and the one genotyped in two individuals with the same
  # genotype); with "remove", every pair in a window is linked
  keep_all <- loci_ld_prune(ld_gt,
    window_size = 10, step_size = 1, thr_r2 = 0.5, min_shared = 61
  )
  no_variance <- apply(show_genotypes(ld_gt), 2, function(g) {
    length(unique(stats::na.omit(g))) < 2
  })
  expect_identical(which(!keep_all), which(no_variance))
  remove_all <- loci_ld_prune(ld_gt,
    window_size = 10, step_size = 1, thr_r2 = 0.5, min_shared = 61,
    insufficient = "remove"
  )
  expect_lt(sum(remove_all), sum(keep_all) / 5)
})

test_that("loci_ld_prune works with several cores and grouped gen_tibbles", {
  one_core <- loci_ld_prune(ld_gt,
    window_size = 10, step_size = 1, thr_r2 = 0.5, n_cores = 1
  )
  expect_identical(
    loci_ld_prune(ld_gt,
      window_size = 10, step_size = 1, thr_r2 = 0.5, n_cores = 2
    ),
    one_core
  )
  grouped_gt <- ld_gt %>%
    mutate(group = rep(c("a", "b"), each = 30)) %>%
    group_by(group)
  expect_identical(
    loci_ld_prune(grouped_gt, window_size = 10, step_size = 1, thr_r2 = 0.5),
    one_core
  )
})
