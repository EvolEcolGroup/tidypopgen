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

test_that("requires a grouped gen_tibble", {
  expect_error(
    test_gt %>% pairwise_pop_fst(method = "Nei87"),
    ".x should be a grouped df"
  )
})

if (rlang::is_installed("hierfstat")) {
  test_that("pairwise_pop_fst Nei87", {
    test_gt <- test_gt %>% dplyr::group_by(population)
    test_hier <- gt_as_hierfstat(test_gt)
    # compare results against hierfstat for Nei87
    nei_gt <- test_gt %>% pairwise_pop_fst(method = "Nei87")

    nei_hier <- hierfstat::pairwise.neifst(test_hier)
    # hiefstat values are rounded to 4 dp
    expect_true(all.equal(
      tidy_dist_matrix(nei_hier)$value,
      round(nei_gt$value, 4)
    ))

    # pair_fst_locus <- test_gt %>% pairwise_pop_fst(by_locus = TRUE) #nolint
  })
}

if (rlang::is_installed("hierfstat")) {
  test_that("pairwise_pop_fst WC84", {
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
      population = c("pop1", "pop1", "pop1", "pop2", "pop2", "pop2", "pop2")
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

    test_gt <- test_gt %>% dplyr::group_by(population)
    wc_tidypopgen <- test_gt %>% pairwise_pop_fst(method = "WC84")
    wc_tidypopgen_per_loc <- test_gt %>%
      pairwise_pop_fst(
        method = "WC84",
        by_locus = TRUE,
        by_locus_type = "matrix"
      )

    # compared to hierfstat
    test_heir <- gt_as_hierfstat(test_gt)
    wc_hierfstat <- hierfstat::pairwise.WCfst(test_heir)
    expect_equal(wc_tidypopgen$value, wc_hierfstat[1, 2])

    # compared to scikit-allel version 1.3.13
    # See create_scikit-allel_test_data for script

    wc_scikit <- as.numeric(readLines(test_path(
      "testdata/fst_scikit-allel",
      "fst_wc.txt"
    )))
    expect_equal(wc_tidypopgen$value, wc_scikit)

    # compare per locus Fst estimate
    wc_scikit_per_loc <- as.numeric(readLines(test_path(
      "testdata/fst_scikit-allel",
      "fst_wc_per_loc.txt"
    )))
    expect_equal(
      wc_scikit_per_loc,
      as.vector(wc_tidypopgen_per_loc$Fst_by_locus)
    )

    ########### test with a monomorphic loci
    test_genotypes <- rbind(
      c(2, 1, 0, 1, 1, 0),
      c(2, 1, 0, NA, 0, 0),
      c(2, NA, 0, 0, 1, 1),
      c(2, 0, 0, 1, 0, 0),
      c(2, 2, 0, 1, 2, 1),
      c(2, 0, 0, 0, NA, 1),
      c(2, 1, 1, 0, 1, NA)
    )
    test_gt <- gen_tibble(
      x = test_genotypes,
      loci = test_loci,
      indiv_meta = test_indiv_meta,
      quiet = TRUE
    )

    test_gt <- test_gt %>% dplyr::group_by(population)
    wc_tidypopgen_mono <- test_gt %>% pairwise_pop_fst(method = "WC84")

    # read in output
    wc_scikit_mono <- as.numeric(readLines(test_path(
      "testdata/fst_scikit-allel",
      "fst_wc_monomorphic.txt"
    )))
    expect_equal(wc_tidypopgen_mono$value, wc_scikit_mono)
    ######################
    # test with 1st loci missing for all individuals in 1st population
    ######################
    test_genotypes <- rbind(
      c(NA, 1, 0, 1, 1, 0),
      c(NA, 1, 0, NA, 0, 0),
      c(NA, NA, 0, 0, 1, 1),
      c(2, 0, 0, 1, 0, 0),
      c(1, 2, 0, 1, 2, 1),
      c(2, 0, 0, 0, NA, 1),
      c(2, 1, 1, 0, 1, NA)
    )
    test_gt <- gen_tibble(
      x = test_genotypes,
      loci = test_loci,
      indiv_meta = test_indiv_meta,
      quiet = TRUE
    )
    test_gt <- test_gt %>% dplyr::group_by(population)
    wc_tidypopgen_fixed_var1 <- test_gt %>% pairwise_pop_fst(method = "WC84")
    # the variant should be ignored
    # so fst should be equal to fst without this variant
    test_genotypes <- rbind(
      c(1, 0, 1, 1, 0),
      c(1, 0, NA, 0, 0),
      c(NA, 0, 0, 1, 1),
      c(0, 0, 1, 0, 0),
      c(2, 0, 1, 2, 1),
      c(0, 0, 0, NA, 1),
      c(1, 1, 0, 1, NA)
    )
    test_loci <- data.frame(
      name = paste0("rs", 1:5),
      chromosome = paste0("chr", c(1, 1, 1, 2, 2)),
      position = as.integer(c(5, 65, 343, 23, 456)),
      genetic_dist = as.double(rep(0, 5)),
      allele_ref = c("T", "C", "G", "C", "T"),
      allele_alt = c("C", NA, "C", "G", "A")
    )
    test_gt <- gen_tibble(
      x = test_genotypes,
      loci = test_loci,
      indiv_meta = test_indiv_meta,
      quiet = TRUE
    )
    test_gt <- test_gt %>% dplyr::group_by(population)
    wc_tidypopgen_missing_var1 <- test_gt %>% pairwise_pop_fst(method = "WC84")
    expect_equal(
      wc_tidypopgen_fixed_var1$value,
      wc_tidypopgen_missing_var1$value
    )
  })
}

test_that("pairwise_pop_fst hudson", {
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
    population = c("pop1", "pop1", "pop1", "pop2", "pop2", "pop2", "pop2")
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

  test_gt <- test_gt %>% dplyr::group_by(population)
  hudson_tidypopgen <- test_gt %>% pairwise_pop_fst(method = "Hudson")
  hudson_tidypopgen_per_loc <- test_gt %>%
    pairwise_pop_fst(
      method = "Hudson", by_locus = TRUE,
      by_locus_type = "matrix"
    )

  # compared to scikit-allel version 1.3.13
  # See create_scikit-allel_test_data for script

  # read in output
  hudson_scikit <- as.numeric(readLines(test_path(
    "testdata/fst_scikit-allel",
    "fst_hudson.txt"
  )))

  # compare total Fst estimate
  hudson_scikit <- as.numeric(readLines(test_path(
    "testdata/fst_scikit-allel",
    "fst_hudson.txt"
  )))
  expect_equal(hudson_tidypopgen$value, hudson_scikit)
  # compare per locus Fst estimate

  hudson_scikit_per_loc <- as.numeric(readLines(test_path(
    "testdata/fst_scikit-allel",
    "fst_hudson_per_loc.txt"
  )))
  expect_equal(
    hudson_scikit_per_loc,
    as.vector(hudson_tidypopgen_per_loc$Fst_by_locus)
  )

  ######################
  # test with a monomorphic loci
  ######################
  test_genotypes <- rbind(
    c(2, 1, 0, 1, 1, 0),
    c(2, 1, 0, NA, 0, 0),
    c(2, NA, 0, 0, 1, 1),
    c(2, 0, 0, 1, 0, 0),
    c(2, 2, 0, 1, 2, 1),
    c(2, 0, 0, 0, NA, 1),
    c(2, 1, 1, 0, 1, NA)
  )
  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  )

  test_gt <- test_gt %>% dplyr::group_by(population)
  hudson_tidypopgen_mono <- test_gt %>% pairwise_pop_fst(method = "Hudson")

  # read in output
  hudson_scikit_mono <- as.numeric(readLines(test_path(
    "testdata/fst_scikit-allel",
    "fst_hudson_monomorphic.txt"
  )))
  expect_equal(hudson_tidypopgen_mono$value, hudson_scikit_mono)

  ######################
  # test with 1st loci missing for all individuals in 1st population
  ######################
  test_genotypes <- rbind(
    c(NA, 1, 0, 1, 1, 0),
    c(NA, 1, 0, NA, 0, 0),
    c(NA, NA, 0, 0, 1, 1),
    c(2, 0, 0, 1, 0, 0),
    c(1, 2, 0, 1, 2, 1),
    c(2, 0, 0, 0, NA, 1),
    c(2, 1, 1, 0, 1, NA)
  )
  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  )
  test_gt <- test_gt %>% dplyr::group_by(population)
  hudson_tidypopgen_fixed_var1 <- test_gt %>%
    pairwise_pop_fst(method = "Hudson")
  # the variant should be ignored
  # so fst should be equal to fst without that variant
  test_genotypes <- rbind(
    c(1, 0, 1, 1, 0),
    c(1, 0, NA, 0, 0),
    c(NA, 0, 0, 1, 1),
    c(0, 0, 1, 0, 0),
    c(2, 0, 1, 2, 1),
    c(0, 0, 0, NA, 1),
    c(1, 1, 0, 1, NA)
  )
  test_loci <- data.frame(
    name = paste0("rs", 1:5),
    chromosome = paste0("chr", c(1, 1, 1, 2, 2)),
    position = as.integer(c(5, 65, 343, 23, 456)),
    genetic_dist = as.double(rep(0, 5)),
    allele_ref = c("T", "C", "G", "C", "T"),
    allele_alt = c("C", NA, "C", "G", "A")
  )
  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  )
  test_gt <- test_gt %>% dplyr::group_by(population)
  hudson_tidypopgen_missing_var1 <- test_gt %>%
    pairwise_pop_fst(method = "Hudson")
  expect_equal(
    hudson_tidypopgen_fixed_var1$value,
    hudson_tidypopgen_missing_var1$value
  )
})

test_that("n_cores can be set", {
  ############
  # Test hudson
  ############
  expect_true(getOption("bigstatsr.check.parallel.blas"))
  test_gt <- test_gt %>% dplyr::group_by(population)
  hudson_one <- test_gt %>% pairwise_pop_fst(method = "Hudson", n_cores = 1)
  hudson_two <- test_gt %>% pairwise_pop_fst(method = "Hudson", n_cores = 2)
  expect_equal(hudson_one$value, hudson_two$value)
  expect_true(getOption("bigstatsr.check.parallel.blas"))

  # test parallel blas is true on exit if function errors
  test_gt <- ungroup(test_gt)
  expect_error(test_gt %>% pairwise_pop_fst(method = "Hudson", n_cores = 2))
  expect_true(getOption("bigstatsr.check.parallel.blas"))

  ############
  # Test nei87
  ############
  expect_true(getOption("bigstatsr.check.parallel.blas"))
  test_gt <- test_gt %>% dplyr::group_by(population)
  nei87_one <- test_gt %>% pairwise_pop_fst(method = "Nei87", n_cores = 1)
  nei87_two <- test_gt %>% pairwise_pop_fst(method = "Nei87", n_cores = 2)
  expect_equal(nei87_one$value, nei87_two$value)
  expect_true(getOption("bigstatsr.check.parallel.blas"))

  # test parallel blas is true on exit if function errors
  test_gt <- ungroup(test_gt)
  expect_error(test_gt %>% pairwise_pop_fst(method = "Nei87", n_cores = 2))
  expect_true(getOption("bigstatsr.check.parallel.blas"))

  ############
  # Test WC84
  ############
  expect_true(getOption("bigstatsr.check.parallel.blas"))
  test_gt <- test_gt %>% dplyr::group_by(population)
  wc84_one <- test_gt %>% pairwise_pop_fst(method = "WC84", n_cores = 1)
  wc84_two <- test_gt %>% pairwise_pop_fst(method = "WC84", n_cores = 2)
  expect_equal(wc84_one$value, wc84_two$value)
  expect_true(getOption("bigstatsr.check.parallel.blas"))

  # test parallel blas is true on exit if function errors
  test_gt <- ungroup(test_gt)
  expect_error(test_gt %>% pairwise_pop_fst(method = "WC84", n_cores = 2))
  expect_true(getOption("bigstatsr.check.parallel.blas"))
})

test_that("type = pairwise", {
  ############
  # Test hudson
  ############
  test_gt <- test_gt %>% dplyr::group_by(population)
  hudson_matrix <-
    test_gt %>% pairwise_pop_fst(method = "Hudson", type = "pairwise")
  hudson_tidy <- test_gt %>% pairwise_pop_fst(method = "Hudson", type = "tidy")
  # check the values are correct
  expect_equal(
    hudson_matrix["pop1", "pop2"],
    subset(
      hudson_tidy,
      hudson_tidy$population_1 == "pop1" &
        hudson_tidy$population_2 == "pop2"
    )$value
  )
  expect_equal(
    hudson_matrix["pop3", "pop1"],
    subset(
      hudson_tidy,
      hudson_tidy$population_1 == "pop1" &
        hudson_tidy$population_2 == "pop3"
    )$value
  )
  ############
  # Test nei87
  ############
  nei87_matrix <-
    test_gt %>% pairwise_pop_fst(method = "Nei87", type = "pairwise")
  nei87_tidy <- test_gt %>% pairwise_pop_fst(method = "Nei87", type = "tidy")
  expect_equal(
    nei87_matrix["pop1", "pop2"],
    subset(
      nei87_tidy,
      nei87_tidy$population_1 == "pop1" &
        nei87_tidy$population_2 == "pop2"
    )$value
  )
  expect_equal(
    nei87_matrix["pop3", "pop1"],
    subset(
      nei87_tidy,
      nei87_tidy$population_1 == "pop1" &
        nei87_tidy$population_2 == "pop3"
    )$value
  )
  ############
  # Test WC84
  ############
  wc84_matrix <-
    test_gt %>% pairwise_pop_fst(method = "WC84", type = "pairwise")
  wc84_tidy <- test_gt %>% pairwise_pop_fst(method = "WC84", type = "tidy")
  expect_equal(
    wc84_matrix["pop1", "pop2"],
    subset(
      wc84_tidy,
      wc84_tidy$population_1 == "pop1" &
        wc84_tidy$population_2 == "pop2"
    )$value
  )
  expect_equal(
    wc84_matrix["pop3", "pop1"],
    subset(
      wc84_tidy,
      wc84_tidy$population_1 == "pop1" &
        wc84_tidy$population_2 == "pop3"
    )$value
  )
})

test_that("by_locus_type", {
  test_gt <- test_gt %>% dplyr::group_by(population)
  ############
  # Test hudson
  ############

  # Matrix
  hudson_locus_matrix <-
    test_gt %>% pairwise_pop_fst(
      method = "Hudson",
      by_locus = TRUE, by_locus_type = "matrix"
    )
  expect_true(is.matrix(hudson_locus_matrix$Fst_by_locus))
  hudson_locus_matrix_pairwise <-
    test_gt %>% pairwise_pop_fst(
      method = "Hudson", by_locus = TRUE,
      type = "pairwise", by_locus_type = "matrix"
    )
  expect_true(is.matrix(hudson_locus_matrix_pairwise$Fst_by_locus))
  expect_equal(
    hudson_locus_matrix$Fst_by_locus,
    hudson_locus_matrix_pairwise$Fst_by_locus
  )
  expect_true(is.matrix(hudson_locus_matrix_pairwise$Fst))

  # List
  hudson_locus_list <-
    test_gt %>% pairwise_pop_fst(
      method = "Hudson",
      by_locus = TRUE, by_locus_type = "list"
    )
  expect_true(is.list(hudson_locus_list$Fst_by_locus))
  hudson_locus_list_pairwise <-
    test_gt %>% pairwise_pop_fst(
      method = "Hudson", by_locus = TRUE,
      type = "pairwise", by_locus_type = "list"
    )
  expect_true(is.list(hudson_locus_list_pairwise$Fst_by_locus))
  expect_equal(
    hudson_locus_list$Fst_by_locus,
    hudson_locus_list_pairwise$Fst_by_locus
  )
  expect_true(is.matrix(hudson_locus_list_pairwise$Fst))

  # Tidy
  hudson_locus_tidy <-
    test_gt %>% pairwise_pop_fst(
      method = "Hudson",
      by_locus = TRUE, by_locus_type = "tidy"
    )
  expect_true(is.data.frame(hudson_locus_tidy$Fst_by_locus))
  hudson_locus_tidy_pairwise <-
    test_gt %>% pairwise_pop_fst(
      method = "Hudson", by_locus = TRUE,
      type = "pairwise", by_locus_type = "tidy"
    )
  expect_true(is.data.frame(hudson_locus_tidy_pairwise$Fst_by_locus))
  expect_equal(
    hudson_locus_tidy$Fst_by_locus,
    hudson_locus_tidy_pairwise$Fst_by_locus
  )
  expect_true(is.matrix(hudson_locus_tidy_pairwise$Fst))

  ############
  # Test nei87
  ############
  # Matrix
  nei87_locus_matrix <-
    test_gt %>% pairwise_pop_fst(
      method = "Nei87",
      by_locus = TRUE, by_locus_type = "matrix"
    )
  expect_true(is.matrix(nei87_locus_matrix$Fst_by_locus))
  nei87_locus_matrix_pairwise <-
    test_gt %>%
    pairwise_pop_fst(
      method = "Nei87", by_locus = TRUE,
      type = "pairwise", by_locus_type = "matrix"
    )
  expect_true(is.matrix(nei87_locus_matrix_pairwise$Fst_by_locus))
  expect_equal(
    nei87_locus_matrix$Fst_by_locus,
    nei87_locus_matrix_pairwise$Fst_by_locus
  )
  expect_true(is.matrix(nei87_locus_matrix_pairwise$Fst))

  # List
  nei_locus_list <-
    test_gt %>% pairwise_pop_fst(
      method = "Nei87",
      by_locus = TRUE, by_locus_type = "list"
    )
  expect_true(is.list(nei_locus_list$Fst_by_locus))
  nei_locus_list_pairwise <-
    test_gt %>% pairwise_pop_fst(
      method = "Nei87", by_locus = TRUE,
      type = "pairwise", by_locus_type = "list"
    )
  expect_true(is.list(nei_locus_list_pairwise$Fst_by_locus))
  expect_equal(
    nei_locus_list$Fst_by_locus,
    nei_locus_list_pairwise$Fst_by_locus
  )
  expect_true(is.matrix(nei_locus_list_pairwise$Fst))

  # Tidy
  nei_locus_tidy <-
    test_gt %>% pairwise_pop_fst(
      method = "Nei87",
      by_locus = TRUE, by_locus_type = "tidy"
    )
  expect_true(is.data.frame(nei_locus_tidy$Fst_by_locus))
  nei_locus_tidy_pairwise <-
    test_gt %>% pairwise_pop_fst(
      method = "Nei87", by_locus = TRUE,
      type = "pairwise", by_locus_type = "tidy"
    )
  expect_true(is.data.frame(nei_locus_tidy_pairwise$Fst_by_locus))
  expect_equal(
    nei_locus_tidy$Fst_by_locus,
    nei_locus_tidy_pairwise$Fst_by_locus
  )
  expect_true(is.matrix(nei_locus_tidy_pairwise$Fst))

  ############
  # Test WC84
  ############
  # Matrix
  wc84_locus_matrix <-
    test_gt %>%
    pairwise_pop_fst(method = "WC84", by_locus = TRUE, by_locus_type = "matrix")
  expect_true(is.matrix(wc84_locus_matrix$Fst_by_locus))
  wc84_locus_matrix_pairwise <-
    test_gt %>%
    pairwise_pop_fst(
      method = "WC84", by_locus = TRUE,
      type = "pairwise", by_locus_type = "matrix"
    )
  expect_true(is.matrix(wc84_locus_matrix_pairwise$Fst_by_locus))
  expect_equal(
    wc84_locus_matrix$Fst_by_locus,
    wc84_locus_matrix_pairwise$Fst_by_locus
  )
  expect_true(is.matrix(wc84_locus_matrix_pairwise$Fst))

  # List
  wc84_locus_list <-
    test_gt %>% pairwise_pop_fst(
      method = "WC84",
      by_locus = TRUE, by_locus_type = "list"
    )
  expect_true(is.list(wc84_locus_list$Fst_by_locus))
  wc84_locus_list_pairwise <-
    test_gt %>% pairwise_pop_fst(
      method = "WC84", by_locus = TRUE,
      type = "pairwise", by_locus_type = "list"
    )
  expect_true(is.list(wc84_locus_list_pairwise$Fst_by_locus))
  expect_equal(
    wc84_locus_list$Fst_by_locus,
    wc84_locus_list_pairwise$Fst_by_locus
  )
  expect_true(is.matrix(wc84_locus_list_pairwise$Fst))

  # Tidy
  wc84_locus_tidy <-
    test_gt %>% pairwise_pop_fst(
      method = "WC84",
      by_locus = TRUE, by_locus_type = "tidy"
    )
  expect_true(is.data.frame(wc84_locus_tidy$Fst_by_locus))
  wc84_locus_tidy_pairwise <-
    test_gt %>% pairwise_pop_fst(
      method = "WC84", by_locus = TRUE,
      type = "pairwise", by_locus_type = "tidy"
    )
  expect_true(is.data.frame(wc84_locus_tidy_pairwise$Fst_by_locus))
  expect_equal(
    wc84_locus_tidy$Fst_by_locus,
    wc84_locus_tidy_pairwise$Fst_by_locus
  )
  expect_true(is.matrix(wc84_locus_tidy_pairwise$Fst))
})
