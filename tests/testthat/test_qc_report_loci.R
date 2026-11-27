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
  # Bonferroni correction, capped at 1
  hwe_res$p_corrected <- apply(
    hwe_res, 1,
    function(row) pmin(1, min(row) * pops)
  )
  # Compare
  expect_equal(loci_report$hwe_p, hwe_res$p_corrected)
  # corrected p-values are valid probabilities
  expect_true(all(loci_report$hwe_p <= 1))
  example_gt <- load_example_gt("gen_tbl") %>% group_by(population)
  example_report <- qc_report_loci(example_gt)
  expect_true(all(example_report$hwe_p >= 0 & example_report$hwe_p <= 1))
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
  # maf (accounting for ploidy) and missingness match the individual functions
  expect_equal(loci_report$maf, loci_maf(test_gt_pseudo))
  expect_equal(loci_report$missingness, loci_missingness(test_gt_pseudo))

  test_gt_pseudo <- test_gt_pseudo %>% group_by(population)

  loci_report_grouped <- qc_report_loci(test_gt_pseudo)
  # grouping does not change maf and missingness
  expect_equal(loci_report_grouped, loci_report)

  # try to plot hwe
  expect_error(autoplot(loci_report, type = "hwe"), "not available")
  expect_error(autoplot(loci_report_grouped, type = "hwe"), "not available")

  # overview plot
  expect_s3_class(autoplot(loci_report, type = "overview"), "patchwork")
  expect_s3_class(autoplot(loci_report_grouped, type = "overview"), "patchwork")

  # all plot
  expect_s3_class(autoplot(loci_report, type = "all"), "ggplot")
  expect_s3_class(autoplot(loci_report_grouped, type = "all"), "ggplot")

  # test other plots individually
})

test_that("qc_report_loci matches loci_maf, loci_missingness and loci_hwe", {
  example_gt <- load_example_gt("gen_tbl")
  loci_report <- suppressMessages(qc_report_loci(example_gt))
  expect_equal(loci_report$snp_id, loci_names(example_gt))
  expect_equal(loci_report$maf, loci_maf(example_gt))
  expect_equal(loci_report$missingness, loci_missingness(example_gt))
  expect_equal(loci_report$hwe_p, loci_hwe(example_gt))
  # codes other than 0, 1, 2 (e.g. dosages) use the individual functions,
  # which must give the same results (bytes 0 to 3 decode as with CODE_012)
  attr(example_gt$genotypes, "fbm")$code256 <- bigsnpr::CODE_DOSAGE
  expect_equal(suppressMessages(qc_report_loci(example_gt)), loci_report)
  grouped_report <- qc_report_loci(example_gt %>% group_by(population))
  attr(example_gt$genotypes, "fbm")$code256 <- bigsnpr::CODE_012
  expect_equal(
    qc_report_loci(example_gt %>% group_by(population)),
    grouped_report
  )
})

test_that("grouped qc_report_loci ignores empty groups", {
  test_gt <- load_example_gt("gen_tbl") %>%
    filter(population == "pop1")
  # a single non-empty group (so that the Bonferroni correction does not cap
  # the p-values at 1), with an empty level before and after it
  test_gt$population <- factor(
    test_gt$population,
    levels = c("pop_empty1", "pop1", "pop_empty2")
  )
  drop_report <- qc_report_loci(test_gt %>% group_by(population))
  keep_gt <- test_gt %>% group_by(population, .drop = FALSE)
  expect_equal(dplyr::n_groups(keep_gt), 3)
  # some loci have a p-value above 0.5, the p-value of the empty groups
  # before they were set to NA
  expect_true(any(drop_report$hwe_p > 0.5))
  expect_equal(qc_report_loci(keep_gt), drop_report)

  # the same when the stats are not computed in a single pass (here with a
  # single individual), where the empty groups of loci_hwe() are NA
  single_gt <- test_gt[1, ]
  drop_report <- qc_report_loci(single_gt %>% group_by(population))
  keep_report <- qc_report_loci(
    single_gt %>% group_by(population, .drop = FALSE)
  )
  expect_false(anyNA(keep_report$hwe_p))
  expect_equal(keep_report, drop_report)
})

test_that("qc_report_loci(hwe = FALSE) matches loci_maf and loci_missingness", {
  example_gt <- load_example_gt("gen_tbl")
  # identical to the two functions, with no message and no hwe_p column
  expect_no_message(loci_report <- qc_report_loci(example_gt, hwe = FALSE))
  expect_s3_class(loci_report, "qc_report_loci")
  expect_named(loci_report, c("snp_id", "maf", "missingness"))
  expect_identical(loci_report$maf, loci_maf(example_gt))
  expect_identical(loci_report$missingness, loci_missingness(example_gt))
  # and to the report with hwe, without its last column
  expect_identical(
    unclass(loci_report),
    unclass(suppressMessages(qc_report_loci(example_gt, n_cores = 1))[1:3])
  )
  # grouping only affects hwe: maf and missingness use all individuals
  expect_identical(
    qc_report_loci(example_gt %>% group_by(population), hwe = FALSE),
    loci_report
  )
  # on a subset of individuals and loci
  sub_gt <- example_gt %>%
    filter(population != "pop1") %>%
    select_loci(c(5, 2, 4))
  sub_report <- qc_report_loci(sub_gt, hwe = FALSE)
  expect_identical(sub_report$maf, loci_maf(sub_gt))
  expect_identical(sub_report$missingness, loci_missingness(sub_gt))
  # dosages are computed separately, with the same results (the FBM is a
  # reference object, so its code is set back afterwards)
  attr(example_gt$genotypes, "fbm")$code256 <- bigsnpr::CODE_DOSAGE
  expect_identical(qc_report_loci(example_gt, hwe = FALSE), loci_report)
  attr(example_gt$genotypes, "fbm")$code256 <- bigsnpr::CODE_012
  # imputed genotypes, used or not
  imputed_gt <- gt_impute_simple(example_gt, method = "mode")
  for (use_imputed in c(TRUE, FALSE)) {
    gt_set_imputed(imputed_gt, use_imputed)
    imputed_report <- qc_report_loci(imputed_gt, hwe = FALSE)
    expect_identical(imputed_report$maf, loci_maf(imputed_gt))
    expect_identical(
      imputed_report$missingness,
      loci_missingness(imputed_gt)
    )
  }
  # (the last report used the observed genotypes, with missing values)
  expect_true(any(imputed_report$missingness > 0))
})
