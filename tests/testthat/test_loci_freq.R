test_that("loci_alt_freq and loci_maf computes correctly", {
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c"),
    population = c("pop1", "pop1", "pop2")
  )
  test_genotypes <- rbind(
    c(1, 1, 0, 1, 1, 2),
    c(2, 1, 0, NA, 0, NA),
    c(2, 2, 0, 0, 1, NA)
  )
  test_loci <- data.frame(
    name = paste0("rs", 1:6),
    chromosome = c(1, 1, 1, 1, 2, 2),
    position = c(3, 5, 65, 343, 23, 456),
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
  # raw frequencies
  freq <- colSums(test_genotypes, na.rm = TRUE) / (c(3, 3, 3, 2, 3, 1) * 2)
  expect_true(all(loci_alt_freq(test_gt$genotypes) == freq))

  # get the counts
  counts <- loci_alt_freq(test_gt, as_counts = TRUE)
  expect_true(all(counts[, 1] / counts[, 2] == freq))

  # convert to minor frequencies
  freq[freq > 0.5] <- 1 - freq[freq > 0.5]
  expect_true(all(loci_maf(test_gt$genotypes) == freq))

  # repeat the tests for a subset of the data
  # remove the 2nd individual and the 3rd and 5th snp
  test_genotypes_subset1 <- test_genotypes[-2, c(-3, -5)]
  test_gt_subset1 <- test_gt %>%
    filter(id != "b") %>%
    select_loci(c(-3, -5))
  freq <- colSums(test_genotypes_subset1, na.rm = TRUE) / (c(2, 2, 2, 1) * 2)
  expect_true(all(loci_alt_freq(test_gt_subset1$genotypes) == freq))

  # convert to minor frequencies
  freq[freq > 0.5] <- 1 - freq[freq > 0.5]
  expect_true(all(loci_maf(test_gt_subset1$genotypes) == freq))

  # repeat the tests for a subset where for loci 6, all genotypes are missing
  # remove the 1st individual and the 3rd and 4th snp
  test_genotypes_subset2 <- test_genotypes[-1, c(-3, -4)]
  test_gt_subset2 <- test_gt %>%
    filter(id != "a") %>%
    select_loci(c(-3, -4))

  # we expect NaN for loci 6 - as both genotypes are NA
  freq <- colSums(test_genotypes_subset2, na.rm = TRUE) / (c(2, 2, 2, NA) * 2)
  expect_equal(loci_alt_freq(test_gt_subset2$genotypes), freq)

  # convert to minor frequencies
  freq[freq > 0.5 & !is.na(freq)] <- 1 - freq[freq > 0.5 & !is.na(freq)]
  expect_equal(loci_maf(test_gt_subset2$genotypes), freq)

  # Test NA in centre
  test_genotypes2 <- rbind(
    c(1, 1, 0, 1, 1, 2),
    c(2, 1, 0, NA, 0, NA),
    c(2, 2, 0, NA, 1, NA)
  )

  test_gt2 <- gen_tibble(
    x = test_genotypes2,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  )

  # repeat the tests for a subset where for loci 6, all genotypes are missing
  # remove the 1st individual
  test_genotypes_subset3 <- test_genotypes2[-1, ]
  test_gt_subset3 <- test_gt2 %>% filter(id != "a") # %>% select_loci(c(-3,-4))

  # we expect NaN for loci 4 and 6 - as both genotypes are NA
  freq2 <- colSums(test_genotypes_subset3, na.rm = TRUE) /
    (c(2, 2, 2, NA, 2, NA) * 2)
  expect_equal(loci_alt_freq(test_gt_subset3$genotypes), freq2)

  # convert to minor frequencies
  freq2[freq2 > 0.5 & !is.na(freq2)] <- 1 - freq2[freq2 > 0.5 & !is.na(freq2)]
  expect_equal(loci_maf(test_gt_subset3$genotypes), freq2)

  # now repeat with multiple blocks of snps
  loci_freq_chunked <- loci_alt_freq(test_gt$genotypes, block_size = 2)
  expect_true(all(loci_alt_freq(test_gt$genotypes) == loci_freq_chunked))
})

test_that("loci_alt_freq and loci_maf on grouped tibbles", {
  test_genotypes <- rbind(
    c(1, 1, 0, 1, 1, 0),
    c(2, 1, 0, NA, 0, 0),
    c(2, NA, 0, 0, 1, 1),
    c(1, 0, 0, 1, 0, 0),
    c(1, 2, 0, 1, 2, 1),
    c(0, 0, 0, 0, NA, 1),
    c(0, 1, 1, 0, 1, NA)
  )
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c", "d", "e", "f", "g"),
    population = c("pop1", "pop1", "pop2", "pop2", "pop1", "pop3", "pop3")
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
  test_gt <- test_gt %>% group_by(population)

  # compute using .grouped_df method
  list <- loci_alt_freq(test_gt, type = "list")
  matrix <- loci_alt_freq(test_gt, type = "matrix")
  tidy <- loci_alt_freq(test_gt, type = "tidy")

  expect_equal(rownames(matrix), show_loci(test_gt)$name)
  expect_equal(colnames(matrix), group_keys(test_gt)$population)
  expect_equal(list[1][[1]], as.vector(matrix[, "pop1"]))
  tidy_pop1 <- tidy %>%
    filter(group == "pop1") %>%
    select(value)
  expect_equal(list[1][[1]], tidy_pop1$value)

  # compute by using group map
  loci_freq_map <- test_gt %>% group_map(.f = ~ loci_alt_freq(.x))
  # use fast cpp code (limit cores to 2)
  loci_freq_grp <- test_gt %>% loci_alt_freq(n_cores = 2)
  all.equal(loci_freq_map, loci_freq_grp)
  # now repeat with multiple blocks of snps
  loci_freq_grp_chunked <- test_gt %>%
    loci_alt_freq(n_cores = 2, block_size = 2)
  expect_true(all.equal(loci_freq_grp, loci_freq_grp_chunked))

  # and now with reframe
  loci_freq_reframe <- test_gt %>% reframe(alt_freq = loci_alt_freq(genotypes))
  loci_freq_direct <- test_gt %>%
    loci_alt_freq() %>%
    arrange(group)
  expect_equal(loci_freq_reframe$alt_freq, loci_freq_direct$value)
  # check that the direct method can take a column genotypes
  loci_freq_direct2 <- test_gt %>%
    loci_alt_freq(genotypes) %>%
    arrange(group)
  expect_equal(loci_freq_reframe$alt_freq, loci_freq_direct2$value)
  # the less intuitive way of reframing
  loci_freq_reframe2 <-
    test_gt %>% reframe(alt_freq = loci_alt_freq(pick(everything())))
  expect_equal(loci_freq_reframe2$alt_freq, loci_freq_reframe$alt_freq)

  # subset
  test_gt_subset <- test_gt %>% select_loci(c(1, 2, 3, 4))
  list <- loci_alt_freq(test_gt_subset, type = "list")
  matrix <- loci_alt_freq(test_gt_subset, type = "matrix")
  tidy <- loci_alt_freq(test_gt_subset, type = "tidy")

  expect_equal(rownames(matrix), show_loci(test_gt_subset)$name)
  expect_equal(colnames(matrix), group_keys(test_gt_subset)$population)
  expect_equal(list[1][[1]], as.vector(matrix[, 1]))
  tidy_pop1 <- tidy %>%
    filter(group == "pop1") %>%
    select(value)
  expect_equal(list[1][[1]], tidy_pop1$value)

  # test a second grouping variable
  test_gt$region <- c("a", "a", "b", "b", "a", "b", "b")
  test_gt <- test_gt %>% group_by(population, region)
  expect_error(
    test_gt %>% loci_alt_freq(),
    "only works with one grouping variable"
  )

  # and now for maf
  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  )
  test_gt <- test_gt %>% group_by(population)

  # compute using .grouped_df method
  list <- loci_maf(test_gt, type = "list")
  matrix <- loci_maf(test_gt, type = "matrix")
  tidy <- loci_maf(test_gt, type = "tidy")

  expect_equal(list[1][[1]], as.vector(matrix[, 1]))
  expect_equal(rownames(matrix), show_loci(test_gt)$name)
  expect_equal(colnames(matrix), group_keys(test_gt)$population)
  tidy_pop1 <- tidy %>%
    filter(group == "pop1") %>%
    select(value)
  expect_equal(list[1][[1]], tidy_pop1$value)

  # compute by using group map
  loci_maf_map <- test_gt %>% group_map(.f = ~ loci_maf(.x))
  # use fast cpp code (limit cores to 2)
  loci_maf_grp <- test_gt %>% loci_maf(n_cores = 2)
  loci_maf_grp_pop1 <- loci_maf_grp %>%
    filter(group == "pop1") %>%
    select(value)
  expect_true(all.equal(loci_maf_map[1][[1]], loci_maf_grp_pop1$value))

  # and now with reframe
  loci_maf_reframe <- test_gt %>% reframe(maf = loci_maf(genotypes))
  loci_maf_direct <- test_gt %>%
    loci_maf() %>%
    arrange(group)
  expect_equal(loci_maf_reframe$maf, loci_maf_direct$value)
  # check that the direct method can take a column genotypes
  loci_maf_direct2 <- test_gt %>%
    loci_maf(genotypes) %>%
    arrange(group)
  expect_equal(loci_maf_reframe$maf, loci_maf_direct2$value)

  # subset
  list <- loci_maf(test_gt_subset, type = "list")
  matrix <- loci_maf(test_gt_subset, type = "matrix")
  tidy <- loci_maf(test_gt_subset, type = "tidy")

  expect_equal(list[1][[1]], as.vector(matrix[, 1]))
  expect_equal(rownames(matrix), show_loci(test_gt_subset)$name)
  expect_equal(colnames(matrix), group_keys(test_gt_subset)$population)
  tidy_pop1 <- tidy %>%
    filter(group == "pop1") %>%
    select(value)
  expect_equal(list[1][[1]], tidy_pop1$value)

  # now repeat with multiple blocks of snps
  loci_freq_grp_chunked <- test_gt %>% loci_maf(n_cores = 2, block_size = 2)
  expect_true(all.equal(loci_maf_grp, loci_freq_grp_chunked))

  # test a second grouping variable
  test_gt$region <- c("a", "a", "b", "b", "a", "b", "b")
  test_gt <- test_gt %>% group_by(population, region)
  expect_error(
    test_gt %>% loci_maf(),
    "only works with one grouping variable"
  )
})
