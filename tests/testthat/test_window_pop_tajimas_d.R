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
test_that("windows_pop_tajimas_d works correctly", {
  window_taj <- windows_pop_tajimas_d(
    type = "list",
    test_gt,
    window_size = 3,
    step_size = 1,
    min_loci = 1
  )
  # subset the tibble to the first population and chr1
  test_gt_chr1 <- test_gt %>%
    select_loci_if(loci_chromosomes(genotypes) == "chr1") %>%
    dplyr::filter(population == "pop1")
  pop_tajimas_d(test_gt_chr1)

  expect_true(
    window_taj[[1]]$stat[1] ==
      pop_tajimas_d(test_gt_chr1)
  )

  # now test 2nd pop and chr2 from position 3 to 5
  test_gt_chr2 <- test_gt %>%
    select_loci_if(
      (loci_chromosomes(genotypes) == "chr2") &
        show_loci(genotypes)$position %in% c(138, 230, 456)
    ) %>%
    dplyr::filter(population == "pop2")
  expect_true(
    window_taj[[2]]$stat[4] ==
      pop_tajimas_d(test_gt_chr2)
  )

  # Additional test for bp-based windows
  test_that("windows_pop_tajimas_d works with bp units", {
    window_taj_bp <- windows_pop_tajimas_d(
      type = "list",
      test_gt,
      window_size = 100,
      step_size = 50,
      size_unit = "bp",
      min_loci = 1
    )

    # The first window should be non-NA
    expect_true(!is.na(window_taj_bp[[1]]$stat[1]))
  })

  # Test for min_loci parameter
  test_that("windows_pop_tajimas_d respects min_loci", {
    # Set min_loci to a value that should produce NA for some windows
    window_taj_high_min <- windows_pop_tajimas_d(
      type = "list",
      test_gt,
      window_size = 2,
      step_size = 1,
      min_loci = 2
    )

    # Check if windows with fewer than min_loci SNPs have NA values
    # TODO a better test would check for specific windows
    expect_true(any(is.na(window_taj_high_min[[1]]$stat)))
  })

  # now ungroupt the gen_tibble
  test_gt <- test_gt %>% ungroup()
  # test the function with no grouping
  window_taj_no_group <- windows_pop_tajimas_d(
    type = "list",
    test_gt,
    window_size = 3,
    step_size = 1,
    min_loci = 1
  )
  expect_true(inherits(window_taj_no_group, "data.frame"))
})

test_that("windows type", {
  window_taj_list <- windows_pop_tajimas_d(
    type = "list",
    test_gt,
    window_size = 3,
    step_size = 1,
    min_loci = 1
  )
  expect_true(is.list(window_taj_list))
  windows_taj_matrix <- windows_pop_tajimas_d(
    type = "matrix",
    test_gt,
    window_size = 3,
    step_size = 1,
    min_loci = 1
  )
  expect_true(inherits(windows_taj_matrix, "data.frame"))
  window_taj_tidy <- windows_pop_tajimas_d(
    type = "tidy",
    test_gt,
    window_size = 3,
    step_size = 1,
    min_loci = 1
  )
  expect_true(inherits(window_taj_tidy, "data.frame"))
  # Compare
  pop1_pop2_tidy <-
    subset(window_taj_tidy, window_taj_tidy$group == "pop3")
  expect_equal(pop1_pop2_tidy$stat, windows_taj_matrix$pop3)
  expect_equal(pop1_pop2_tidy$stat, window_taj_list$pop3$stat)
  expect_equal(windows_taj_matrix$pop3, window_taj_list$pop3$stat)
})
