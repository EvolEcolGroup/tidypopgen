test_that("grouped hwe within qc_report_loci", {
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
      "pop1", "pop1", "pop2", "pop2", "pop1", "pop2", "pop1",
      "pop2", "pop2", "pop1"
    )
  )
  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  )
  test_gt <- test_gt %>% dplyr::group_by(population)

  loci_report <- qc_report_loci(test_gt)

  # Calculate with group_map
  hwe_res <- test_gt %>% group_map(.f = ~ loci_hwe(.x$genotypes))
  hwe_res <- do.call("cbind", hwe_res)
  hwe_res <- as.data.frame(hwe_res)
  pops <- test_gt %>%
    select(dplyr::group_vars(test_gt)) %>%
    dplyr::pull(1)
  pops <- length(unique(pops))
  hwe_res$p_corrected <- apply(hwe_res, 1, function(row) min(row) * pops)
  # Compare
  expect_equal(loci_report$hwe_p, hwe_res$p_corrected)
})
