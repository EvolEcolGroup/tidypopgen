test_genotypes <- rbind(
  c(1, 1, 0, 1, 1, 0),
  c(2, 2, 0, NA, 0, 0),
  c(2, NA, 0, 0, 1, 1),
  c(2, 0, 0, 2, 0, 0),
  c(1, 2, 0, 1, 2, 1),
  c(0, 0, 0, 0, NA, 2),
  c(0, 1, 2, 0, 1, NA)
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

test_that("gt_pseudohaploid correctly deals with ploidy", {
  # confirm that the gen_tibble is diploid for the moment
  expect_equal(show_ploidy(test_gt), 2)
  # now use gt_pseudohaploid
  test_gt_pseudo <- gt_pseudohaploid(test_gt)
  # now check that ploidy is -2
  expect_equal(show_ploidy(test_gt_pseudo), -2)
  ## check the individual ploidies
  expect_equal(
    indiv_ploidy(test_gt_pseudo),
    c(2L, 1L, 2L, 1L, 2L, 1L, 2L)
  )

  # check that we can reprocess it (use a different number of loci to get
  # different result)
  test_gt_pseudo2 <- gt_pseudohaploid(test_gt_pseudo, test_n_loci = 4)
  expect_equal(
    indiv_ploidy(test_gt_pseudo2),
    c(2L, 1L, 1L, 1L, 2L, 1L, 2L)
  )

  ## now test which functions work and which fail with pseudohaploid data (we
  # should get the right frequencies (grouped and ungrouped, and pairwise pop
  # fst)) missingness should also work but most other indiv stats will fail, and
  # so will pop estimates that requires heterozygote counts
  expect_error(
    indiv_het_obs(test_gt_pseudo),
    "this function only works on diploid data"
  )
  expect_error(
    indiv_inbreeding(test_gt_pseudo),
    "this function only works on diploid data"
  )
  expect_error(
    loci_hwe(test_gt_pseudo),
    "this function only works on diploid data"
  )
  expect_error(
    pop_fis(test_gt_pseudo),
    "this function only works on diploid data"
  )
  expect_error(
    pop_fst(test_gt_pseudo),
    "this function only works on diploid data"
  )
  # now check frequencies
  pseudo_freq <- loci_alt_freq(test_gt_pseudo)
  ploidy_hap <- c(2L, 1L, 2L, 1L, 2L, 1L, 2L)
  # convert genotypes to alt counts by ploidy
  test_genotypes_hap <- sweep(test_genotypes, 1,
    3 - ploidy_hap,
    FUN = "/"
  )
  test_valid_alleles <- matrix(ploidy_hap,
    ncol = ncol(test_genotypes_hap),
    nrow = nrow(test_genotypes_hap)
  )
  test_valid_alleles[is.na(test_genotypes_hap)] <- NA
  test_valid_alleles <- colSums(test_valid_alleles, na.rm = TRUE)
  expect_identical(
    pseudo_freq,
    colSums(test_genotypes_hap, na.rm = TRUE) /
      test_valid_alleles
  )
  counts <- loci_alt_freq(test_gt_pseudo, as_counts = TRUE)
  expect_true(
    all(
      cbind(
        colSums(test_genotypes_hap, na.rm = TRUE),
        test_valid_alleles
      ) ==
        loci_alt_freq(test_gt_pseudo, as_counts = TRUE)
    )
  )

  # and missingness
  expect_true(
    all(
      loci_missingness(test_gt_pseudo) ==
        colSums(is.na(test_genotypes_hap)) / nrow(test_genotypes_hap)
    )
  )

  # and maf
  pseudo_maf <- loci_maf(test_gt_pseudo)
  pseudo_freq[pseudo_freq > 0.5] <- 1 - pseudo_freq[pseudo_freq > 0.5]
  expect_identical(pseudo_maf, pseudo_freq)
})



test_that("gt_pseudohaploid on grouped tibble correctly deals with ploidy", {
  ## now group it
  test_gt <- test_gt %>% group_by(population)
  # confirm that the gen_tibble is diploid for the moment
  expect_equal(show_ploidy(test_gt), 2)
  # now use gt_pseudohaploid
  test_gt_pseudo <- gt_pseudohaploid(test_gt)
  # now check that ploidy is -2
  expect_equal(show_ploidy(test_gt_pseudo), -2)
  ## check the individual ploidies
  expect_equal(
    indiv_ploidy(test_gt_pseudo),
    c(2L, 1L, 2L, 1L, 2L, 1L, 2L)
  )

  # check that we can reprocess it (use a different number of loci to get
  # different result)
  test_gt_pseudo2 <- gt_pseudohaploid(test_gt_pseudo, test_n_loci = 4)
  expect_equal(
    indiv_ploidy(test_gt_pseudo2),
    c(2L, 1L, 1L, 1L, 2L, 1L, 2L)
  )

  ## now test which functions work and which fail with pseudohaploid data (we
  # should get the right frequencies (grouped and ungrouped, and pairwise pop
  # fst)) missingness should also work but most other indiv stats will fail, and
  # so will pop estimates that requires heterozygote counts
  expect_error(
    indiv_het_obs(test_gt_pseudo),
    "this function only works on diploid data"
  )
  expect_error(
    indiv_inbreeding(test_gt_pseudo),
    "this function only works on diploid data"
  )
  expect_error(
    loci_hwe(test_gt_pseudo),
    "this function only works on diploid data"
  )
  expect_error(
    pop_fis(test_gt_pseudo),
    "this function only works on diploid data"
  )
  expect_error(
    pop_fst(test_gt_pseudo),
    "this function only works on diploid data"
  )
  # now check frequencies with reframe
  loci_freq_reframe <- test_gt %>% reframe(alt_freq = loci_alt_freq(genotypes))
  loci_freq_direct <- test_gt %>%
    loci_alt_freq() %>%
    arrange(group)
  expect_equal(loci_freq_reframe$alt_freq, loci_freq_direct$value)

  # maf
  pseudo_maf <- test_gt %>% reframe(maf = loci_maf(genotypes))
  pseudo_maf_direct <- test_gt %>%
    loci_maf() %>%
    arrange(group)
  expect_equal(pseudo_maf$maf, pseudo_maf_direct$value)

  # apply Fst on the pseudohap gen_tibble
  pseudo_fst <- pairwise_pop_fst(test_gt_pseudo)
  # confirm the results change
  expect_false(all(pseudo_fst$value == pairwise_pop_fst(test_gt)$value))
  # TODO we should test the values are correct

  # only Hudson works on pseudohaploids
  expect_error(pairwise_pop_fst(test_gt_pseudo, method = "Nei87"),
    "only `method = Hudson` is valid",
    fixed = TRUE
  )

  # PBS also works if we use Hudson
  pseudo_pbs <- nwise_pop_pbs(test_gt_pseudo)
  expect_false(identical(
    pseudo_pbs,
    nwise_pop_pbs(test_gt)
  ))
  # only Hudson works on pseudohaploids
  expect_error(nwise_pop_pbs(test_gt_pseudo, fst_method = "Nei87"),
    "only `method = Hudson` is valid",
    fixed = TRUE
  )
})

test_that("we can rbind pseudohaploids and diploids", {
  # create a small diploid gen_tibble
  test_genotypes_dip <- rbind(
    c(1, 1, 0, 1, NA, 0),
    c(2, 1, 0, 0, 0, 0),
    c(1, NA, 0, 0, 2, 1)
  )
  test_indiv_meta_dip <- data.frame(
    id = c("m", "n", "o"),
    population = c("pop4", "pop4", "pop5")
  )

  test_gt_dip <- gen_tibble(
    x = test_genotypes_dip,
    loci = test_loci,
    indiv_meta = test_indiv_meta_dip,
    quiet = TRUE
  )

  test_gt_pseudo <- gt_pseudohaploid(test_gt)

  merged_gt <- rbind(test_gt_pseudo, test_gt_dip,
    quiet = TRUE,
    backingfile = tempfile()
  )
  expect_true(show_ploidy(merged_gt) == -2)
  expect_true(all(indiv_ploidy(merged_gt) == c(
    indiv_ploidy(test_gt_pseudo),
    indiv_ploidy(test_gt_dip)
  )))
})

test_that("gt_pseudohaploid updates ploidy after dropping pseudohploids", {
  test_genotypes <- rbind(
    c(2, NA, 0, NA, 2, 0),
    c(2, 2, 0, NA, 0, NA),
    c(2, NA, 0, 0, 1, 1),
    c(2, 0, 0, 1, 0, 0),
    c(1, 2, 0, 1, 2, 1),
    c(0, 0, 1, 0, 1, 2),
    c(0, 1, 2, 0, 1, 1)
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

  # before gt_pseudohaploid ploidy is 2
  expect_equal(test_gt %>% show_ploidy(), 2)
  # now use gt_pseudohaploid
  test_gt <- gt_pseudohaploid(test_gt)
  # now check that ploidy is -2
  expect_equal(test_gt %>% show_ploidy(), -2)
  # and functions requiring diploid data throw errors
  expect_error(
    round(indiv_het_obs(test_gt), 1),
    "this function only works on diploid data"
  )

  # Now try filtering out pseudohaploid individuals
  test_gt <- test_gt %>% filter(indiv_missingness(genotypes) < 0.15)

  # Functions using stopifnot_diploid should now work
  # becasue data no longer contain pseudohaploids
  expect_equal(length(indiv_het_obs(test_gt)), 4)
  expect_equal(round(indiv_het_obs(test_gt), 1), c(0.2, 0.5, 0.3, 0.5))

  # Rerun gt_pseudohaploid to update gen_tibble
  test_gt <- gt_pseudohaploid(test_gt)
  # check individuals
  expect_equal(test_gt %>% indiv_ploidy(), rep(2, 4))
  # check ploidy
  expect_equal(test_gt %>% show_ploidy(), 2)
})
