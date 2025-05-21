test_genotypes <- rbind(
  c(1, 0, 1, 1, 0, 2, 1, 0),
  c(1, 0, NA, 0, 0, 0, 1, 2),
  c(NA, 0, 0, 1, 1, 0, 1, 0),
  c(0, 0, 1, 0, 0, 0, 1, 0),
  c(2, 0, 1, 2, 1, 1, 2, 2),
  c(0, 0, 0, NA, 1, 0, 0, 1),
  c(1, 1, 0, 1, NA, 0, 1, 0),
  c(0, 0, 1, 0, 0, 0, 1, 2),
  c(1, 0, 1, 0, 0, 0, 1, 0),
  c(0, 0, 1, 0, 0, 0, 1, 2)
)
test_loci <- data.frame(
  name = paste0("rs", 1:8),
  chromosome = paste0("chr", c(1, 1, 1, 2, 2, 2, 2, 2)),
  position = as.integer(c(5, 65, 343, 23, 56, 138, 230, 456)),
  genetic_dist = as.double(rep(0, 8)),
  allele_ref = c("T", "C", "G", "C", "T", "A", "C", "G"),
  allele_alt = c("C", NA, "C", "G", "A", "T", "A", "C")
)
test_indiv_meta <- data.frame(
  id = c("a", "b", "c", "d", "e", "f", "g", "h", "i", "l"),
  population = c(
    "pop1", "pop1", "pop2", "pop2", "pop1", "pop3", "pop3",
    "pop2", "pop3", "pop3"
  )
)
test_gt <- gen_tibble(
  x = test_genotypes,
  loci = test_loci,
  indiv_meta = test_indiv_meta,
  quiet = TRUE
)
test_gt <- test_gt %>% dplyr::group_by(population)

# testing the infrastructure for windowing
test_that("pairwise_pop_fst num and dem are returned correctly", {
  hudson_gt_fst <- test_gt %>%
    pairwise_pop_fst(
      method = "Hudson",
      by_locus = TRUE,
      by_locus_type = "matrix"
    )
  expect_message(
    hudson_gt_fst_nd <- test_gt %>%
      pairwise_pop_fst(
        method = "Hudson",
        return_num_dem = TRUE
      ),
    "`by_locus` set to TRUE because `return_num_dem = TRUE`"
  )
  # check that the nums and denominators are indeed the right ones to
  # generate the Fst for each locus
  hudson_gt_from_num_dem <- hudson_gt_fst_nd$Fst_by_locus_num /
    hudson_gt_fst_nd$Fst_by_locus_den
  colnames(hudson_gt_from_num_dem) <- paste0(
    "fst_",
    colnames(hudson_gt_from_num_dem)
  )
  expect_equal(hudson_gt_from_num_dem, hudson_gt_fst$Fst_by_locus)

  # repeat the same for WC84
  wc_gt_fst <- test_gt %>%
    pairwise_pop_fst(
      method = "WC84",
      by_locus = TRUE,
      by_locus_type = "matrix"
    )
  expect_message(
    wc_gt_fst_nd <- test_gt %>%
      pairwise_pop_fst(
        method = "WC84",
        return_num_dem = TRUE
      ),
    "`by_locus` set to TRUE because `return_num_dem = TRUE`"
  )
  # check that the nums and denominators are indeed the right ones to
  # generate the Fst for each locus
  wc_gt_from_num_dem <- wc_gt_fst_nd$Fst_by_locus_num /
    wc_gt_fst_nd$Fst_by_locus_den
  colnames(wc_gt_from_num_dem) <- paste0(
    "fst_",
    colnames(wc_gt_from_num_dem)
  )
  expect_equal(wc_gt_from_num_dem, wc_gt_fst$Fst_by_locus)

  # And for Nei
  nei_gt_fst <- test_gt %>%
    pairwise_pop_fst(
      method = "Nei87",
      by_locus = TRUE,
      by_locus_type = "matrix"
    )
  expect_message(
    nei_gt_fst_nd <- test_gt %>%
      pairwise_pop_fst(
        method = "Nei87",
        return_num_dem = TRUE
      ),
    "`by_locus` set to TRUE because `return_num_dem = TRUE`"
  )
  # check that the nums and denominators are indeed the right ones to
  # generate the Fst for each locus
  nei_gt_from_num_dem <- nei_gt_fst_nd$Fst_by_locus_num /
    nei_gt_fst_nd$Fst_by_locus_den
  colnames(nei_gt_from_num_dem) <- paste0(
    "fst_",
    colnames(nei_gt_from_num_dem)
  )
  expect_equal(nei_gt_from_num_dem, nei_gt_fst$Fst_by_locus)
})

test_that("windows_pairwise_pop_fst works correctly", {
  hudson_gt_fst_nd <- test_gt %>%
    pairwise_pop_fst(
      method = "Hudson",
      return_num_dem = TRUE,
      by_locus = TRUE
    )
  # check that pop Fst by SNP is calculated correctly
  snp_window <- test_gt %>%
    windows_pairwise_pop_fst(
      window_size = 3,
      step_size = 2,
      size_unit = "snp",
      min_loci = 2
    )
  # check step size
  expect_equal(
    snp_window$start[3] - snp_window$start[2],
    2
  )
  # check chrm 1, window 1 Fst calculation
  # find first 3 SNP loci
  ch1_wind1_pos <- test_gt %>%
    show_loci() %>%
    filter(chromosome == "chr1") %>%
    slice_head(n = 3)
  # pairwise fst for first 3 SNPs
  ch1_wind1_fst <- test_gt %>%
    select_loci(ch1_wind1_pos$big_index) %>%
    pairwise_pop_fst(
      method = "Hudson",
      by_locus = TRUE
    )
  # pop1 vs pop2
  expect_equal(
    snp_window$fst_pop1.pop2[1],
    ch1_wind1_fst$Fst$value[1]
  )
  # pop1 vs pop3
  expect_equal(
    snp_window$fst_pop1.pop3[1],
    ch1_wind1_fst$Fst$value[2]
  )
  # pop2 vs pop3
  expect_equal(
    snp_window$fst_pop2.pop3[1],
    ch1_wind1_fst$Fst$value[3]
  )
  # manually calculate chrm 1, window 1 Fst
  manual_pop1_pop2 <- mean(hudson_gt_fst_nd$Fst_by_locus_num[1:3, 1]) /
    mean(hudson_gt_fst_nd$Fst_by_locus_den[1:3, 1])
  manual_pop1_pop3 <- mean(hudson_gt_fst_nd$Fst_by_locus_num[1:3, 2]) /
    mean(hudson_gt_fst_nd$Fst_by_locus_den[1:3, 2])
  manual_pop2_pop3 <- mean(hudson_gt_fst_nd$Fst_by_locus_num[1:3, 3]) /
    mean(hudson_gt_fst_nd$Fst_by_locus_den[1:3, 3])
  # pop1 vs pop2
  expect_equal(
    snp_window$fst_pop1.pop2[1],
    manual_pop1_pop2
  )
  # pop1 vs pop3
  expect_equal(
    snp_window$fst_pop1.pop3[1],
    manual_pop1_pop3
  )
  # pop2 vs pop3
  expect_equal(
    snp_window$fst_pop2.pop3[1],
    manual_pop2_pop3
  )
  # check chrm 2, window 1 Fst calculation
  # find first 3 SNP loci on chr 2
  ch2_wind1_pos <- test_gt %>%
    show_loci() %>%
    filter(chromosome == "chr2") %>%
    slice_head(n = 3)
  # pairwise fst for first 3 SNPs chr 2
  ch2_wind1_fst <- test_gt %>%
    select_loci(ch2_wind1_pos$big_index) %>%
    pairwise_pop_fst(
      method = "Hudson",
      by_locus = TRUE
    )
  # pop1 vs pop2
  expect_equal(
    snp_window$fst_pop1.pop2[2],
    ch2_wind1_fst$Fst$value[1]
  )
  # pop1 vs pop3
  expect_equal(
    snp_window$fst_pop1.pop3[2],
    ch2_wind1_fst$Fst$value[2]
  )
  # pop2 vs pop3
  expect_equal(
    snp_window$fst_pop2.pop3[2],
    ch2_wind1_fst$Fst$value[3]
  )
  # manually calculate chrm 2, window 1 Fst
  manual_pop1_pop2 <- mean(hudson_gt_fst_nd$Fst_by_locus_num[4:6, 1]) /
    mean(hudson_gt_fst_nd$Fst_by_locus_den[4:6, 1])
  manual_pop1_pop3 <- mean(hudson_gt_fst_nd$Fst_by_locus_num[4:6, 2]) /
    mean(hudson_gt_fst_nd$Fst_by_locus_den[4:6, 2])
  manual_pop2_pop3 <- mean(hudson_gt_fst_nd$Fst_by_locus_num[4:6, 3]) /
    mean(hudson_gt_fst_nd$Fst_by_locus_den[4:6, 3])
  # pop1 vs pop2
  expect_equal(
    snp_window$fst_pop1.pop2[2],
    manual_pop1_pop2
  )
  # pop1 vs pop3
  expect_equal(
    snp_window$fst_pop1.pop3[2],
    manual_pop1_pop3
  )
  # pop2 vs pop3
  expect_equal(
    snp_window$fst_pop2.pop3[2],
    manual_pop2_pop3
  )
  # check chrm 2, window 2 Fst calculation
  # find 3-5 SNP loci on chr 2
  ch2_wind2_pos <- test_gt %>%
    show_loci() %>%
    filter(chromosome == "chr2") %>%
    slice_tail(n = 3)
  # pairwise fst for 3-5 SNPs chr 2
  ch2_wind2_fst <- test_gt %>%
    select_loci(ch2_wind2_pos$big_index) %>%
    pairwise_pop_fst(
      method = "Hudson",
      by_locus = TRUE
    )
  # pop1 vs pop2
  expect_equal(
    snp_window$fst_pop1.pop2[3],
    ch2_wind2_fst$Fst$value[1]
  )
  # pop1 vs pop3
  expect_equal(
    snp_window$fst_pop1.pop3[3],
    ch2_wind2_fst$Fst$value[2]
  )
  # pop2 vs pop3
  expect_equal(
    snp_window$fst_pop2.pop3[3],
    ch2_wind2_fst$Fst$value[3]
  )
  # manually calculate chrm 2, window 2 Fst
  manual_pop1_pop2 <- mean(hudson_gt_fst_nd$Fst_by_locus_num[6:8, 1]) /
    mean(hudson_gt_fst_nd$Fst_by_locus_den[6:8, 1])
  manual_pop1_pop3 <- mean(hudson_gt_fst_nd$Fst_by_locus_num[6:8, 2]) /
    mean(hudson_gt_fst_nd$Fst_by_locus_den[6:8, 2])
  manual_pop2_pop3 <- mean(hudson_gt_fst_nd$Fst_by_locus_num[6:8, 3]) /
    mean(hudson_gt_fst_nd$Fst_by_locus_den[6:8, 3])
  # pop1 vs pop2
  expect_equal(
    snp_window$fst_pop1.pop2[3],
    manual_pop1_pop2
  )
  # pop1 vs pop3
  expect_equal(
    snp_window$fst_pop1.pop3[3],
    manual_pop1_pop3
  )
  # pop2 vs pop3
  expect_equal(
    snp_window$fst_pop2.pop3[3],
    manual_pop2_pop3
  )
  # check that pop Fst by BP is calculated correctly
  bp_window <- windows_pairwise_pop_fst(test_gt,
    window_size = 200,
    step_size = 100,
    size_unit = "bp",
    min_loci = 1
  )
  # check step size
  expect_equal(
    bp_window$start[2] - bp_window$start[1],
    100
  )
  # check chrm 1, window 1 Fst calculation
  # find loci between positions 1-200
  ch1_wind1_pos <- test_gt %>%
    show_loci() %>%
    filter(chromosome == "chr1" & position >= 1 & position <= 200) %>%
    select(big_index)
  # pairwise fst for 1-200 positions
  ch1_wind1_fst <- test_gt %>%
    select_loci(ch1_wind1_pos$big_index) %>%
    pairwise_pop_fst(
      method = "Hudson",
      by_locus = TRUE
    )
  # pop1 vs pop2
  expect_equal(
    bp_window$fst_pop1.pop2[1],
    ch1_wind1_fst$Fst$value[1]
  )
  # pop1 vs pop3
  expect_equal(
    bp_window$fst_pop1.pop3[1],
    ch1_wind1_fst$Fst$value[2]
  )
  # pop2 vs pop3
  expect_equal(
    bp_window$fst_pop2.pop3[1],
    ch1_wind1_fst$Fst$value[3]
  )
  # manually calculate chrm 1, window 1 Fst
  manual_pop1_pop2 <- mean(hudson_gt_fst_nd$Fst_by_locus_num[1:2, 1]) /
    mean(hudson_gt_fst_nd$Fst_by_locus_den[1:2, 1])
  manual_pop1_pop3 <- mean(hudson_gt_fst_nd$Fst_by_locus_num[1:2, 2]) /
    mean(hudson_gt_fst_nd$Fst_by_locus_den[1:2, 2])
  manual_pop2_pop3 <- mean(hudson_gt_fst_nd$Fst_by_locus_num[1:2, 3]) /
    mean(hudson_gt_fst_nd$Fst_by_locus_den[1:2, 3])
  # pop1 vs pop2
  expect_equal(
    bp_window$fst_pop1.pop2[1],
    manual_pop1_pop2
  )
  # pop1 vs pop3
  expect_equal(
    bp_window$fst_pop1.pop3[1],
    manual_pop1_pop3
  )
  # pop2 vs pop3
  expect_equal(
    bp_window$fst_pop2.pop3[1],
    manual_pop2_pop3
  )
  # check chrm 1, window 2 (min_loci < 1) is NA
  expect_true(is.na(bp_window$fst_pop1.pop2[2]))
  expect_true(is.na(bp_window$fst_pop1.pop3[2]))
  expect_true(is.na(bp_window$fst_pop2.pop3[2]))

  # check chrm 2, window 1 Fst calculation
  # find loci between positions 1-200 on chr 2
  ch2_wind1_pos <- test_gt %>%
    show_loci() %>%
    filter(chromosome == "chr2" & position >= 1 & position <= 200) %>%
    select(big_index)
  # pairwise fst for 1-200 positions chr 2
  ch2_wind1_fst <- test_gt %>%
    select_loci(ch2_wind1_pos$big_index) %>%
    pairwise_pop_fst(
      method = "Hudson",
      by_locus = TRUE
    )
  # pop1 vs pop2
  expect_equal(
    bp_window$fst_pop1.pop2[4],
    ch2_wind1_fst$Fst$value[1]
  )
  # pop1 vs pop3
  expect_equal(
    bp_window$fst_pop1.pop3[4],
    ch2_wind1_fst$Fst$value[2]
  )
  # pop2 vs pop3
  expect_equal(
    bp_window$fst_pop2.pop3[4],
    ch2_wind1_fst$Fst$value[3]
  )
  # manually calculate chrm 2, window 1 Fst
  manual_pop1_pop2 <- mean(hudson_gt_fst_nd$Fst_by_locus_num[4:6, 1]) /
    mean(hudson_gt_fst_nd$Fst_by_locus_den[4:6, 1])
  manual_pop1_pop3 <- mean(hudson_gt_fst_nd$Fst_by_locus_num[4:6, 2]) /
    mean(hudson_gt_fst_nd$Fst_by_locus_den[4:6, 2])
  manual_pop2_pop3 <- mean(hudson_gt_fst_nd$Fst_by_locus_num[4:6, 3]) /
    mean(hudson_gt_fst_nd$Fst_by_locus_den[4:6, 3])
  # pop1 vs pop2
  expect_equal(
    bp_window$fst_pop1.pop2[4],
    manual_pop1_pop2
  )
  # pop1 vs pop3
  expect_equal(
    bp_window$fst_pop1.pop3[4],
    manual_pop1_pop3
  )
  # pop2 vs pop3
  expect_equal(
    bp_window$fst_pop2.pop3[4],
    manual_pop2_pop3
  )
  # check chrm 2, window 5 Fst calculation
  # find loci between positions 401-600 on chr 2
  ch2_wind5_pos <- test_gt %>%
    show_loci() %>%
    filter(chromosome == "chr2" & position >= 401 & position <= 600) %>%
    select(big_index)
  # pairwise fst for 401-600 positions chr 2
  ch2_wind5_fst <- test_gt %>%
    select_loci(ch2_wind5_pos$big_index) %>%
    pairwise_pop_fst(
      method = "Hudson",
      by_locus = TRUE
    )
  # pop1 vs pop2
  expect_equal(
    bp_window$fst_pop1.pop2[8],
    ch2_wind5_fst$Fst$value[1]
  )
  # pop1 vs pop3
  expect_equal(
    bp_window$fst_pop1.pop3[8],
    ch2_wind5_fst$Fst$value[2]
  )
  # pop2 vs pop3
  expect_equal(
    bp_window$fst_pop2.pop3[8],
    ch2_wind5_fst$Fst$value[3]
  )
  # manually calculate chrm 2, window 5 Fst
  manual_pop1_pop2 <- mean(hudson_gt_fst_nd$Fst_by_locus_num[8, 1]) /
    mean(hudson_gt_fst_nd$Fst_by_locus_den[8, 1])
  manual_pop1_pop3 <- mean(hudson_gt_fst_nd$Fst_by_locus_num[8, 2]) /
    mean(hudson_gt_fst_nd$Fst_by_locus_den[8, 2])
  manual_pop2_pop3 <- mean(hudson_gt_fst_nd$Fst_by_locus_num[8, 3]) /
    mean(hudson_gt_fst_nd$Fst_by_locus_den[8, 3])
  # pop1 vs pop2
  expect_equal(
    bp_window$fst_pop1.pop2[8],
    manual_pop1_pop2
  )
  # pop1 vs pop3
  expect_equal(
    bp_window$fst_pop1.pop3[8],
    manual_pop1_pop3
  )
  # pop2 vs pop3
  expect_equal(
    bp_window$fst_pop2.pop3[8],
    manual_pop2_pop3
  )
})

test_that("windows type", {
  fst_windows_matrix <- test_gt %>%
    windows_pairwise_pop_fst(
      type = "matrix",
      window_size = 3,
      step_size = 2,
      size_unit = "snp",
      min_loci = 2
    )
  expect_true(is.data.frame(fst_windows_matrix))
  fst_windows_tidy <- test_gt %>%
    windows_pairwise_pop_fst(
      type = "tidy",
      window_size = 3,
      step_size = 2,
      size_unit = "snp",
      min_loci = 2
    )
  expect_true(is.data.frame(fst_windows_tidy))
  # Compare
  pop1_pop2_tidy <-
    subset(fst_windows_tidy, fst_windows_tidy$stat_name == "fst_pop1.pop2")
  expect_equal(pop1_pop2_tidy$value, fst_windows_matrix$fst_pop1.pop2)
})
