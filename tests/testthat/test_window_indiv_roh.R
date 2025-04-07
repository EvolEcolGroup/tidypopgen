skip_if(!requireNamespace("detectRUNS", quietly = TRUE))

test_that("window_indiv_roh works without backingfile update", {
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c"),
    population = c("pop1", "pop1", "pop2")
  )
  test_genotypes <- rbind(
    c(1, 1, 1, 1, 1, 1, 0),
    c(2, 1, 0, 0, 0, 0, 1),
    c(2, 2, 0, 0, 1, 1, 2)
  )
  test_loci <- data.frame(
    name = paste0("rs", 1:7),
    chromosome = paste0("chr", c(1, 1, 1, 1, 2, 2, 2)),
    position = as.integer(c(1000, 1500, 2000, 2500, 50000, 60000, 60500)),
    genetic_dist = as.double(rep(0, 7)),
    allele_ref = c("A", "T", "C", "G", "C", "T", "T"),
    allele_alt = c("T", "C", NA, "C", "G", "A", "A")
  )

  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  )

  suppressMessages(test_roh1 <- window_indiv_roh(
    test_gt,
    window_size = 4,
    min_snp = 2,
    min_density = 1 / 500,
    max_gap = 5000,
    min_length_bps = 1000,
    threshold = 0.05
  ))

  # remove one snp
  test_gt <- test_gt %>% select_loci(-7)

  # Check roh still works
  suppressMessages(test_roh2 <- window_indiv_roh(
    test_gt,
    window_size = 4,
    min_snp = 2,
    min_density = 1 / 500,
    max_gap = 5000,
    min_length_bps = 1000,
    threshold = 0.05
  ))

  expect_equal(test_roh1, test_roh2)
})
