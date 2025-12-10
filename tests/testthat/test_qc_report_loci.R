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

test_that("qc_report_loci for pseudohaploid data", {
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

  test_gt_pseudo <- gt_pseudohaploid(test_gt)

  loci_report <- qc_report_loci(test_gt_pseudo)

  test_gt_pseudo <- test_gt_pseudo %>% group_by(population)

  loci_report_grouped <- qc_report_loci(test_gt_pseudo)

  # try to plot hwe
  expect_error(autoplot(loci_report, type = "hwe"), "not available")
  expect_error(autoplot(loci_report_grouped, type = "hwe"), "not available")

  # overview plot
  expect_s3_class(autoplot(loci_report, type = "overview"), "upset")
  expect_s3_class(autoplot(loci_report_grouped, type = "overview"), "upset")

  # all plot
  expect_s3_class(autoplot(loci_report, type = "all"), "ggplot")
  expect_s3_class(autoplot(loci_report_grouped, type = "all"), "ggplot")

  # test other plots individually
})
