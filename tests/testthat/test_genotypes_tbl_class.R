test_indiv_meta <- data.frame(
  id = c("a", "b", "c", "d", "e", "f"),
  population = c("pop1", "pop1", "pop1", "pop2", "pop2", "pop2")
)
test_genotypes <- rbind(
  c(1, 1, 0, 1, 1, 0),
  c(2, 1, 0, 0, 0, 0),
  c(2, 2, 0, 0, 1, 1),
  c(0, 1, 2, 1, 0, 0),
  c(1, 0, 2, 2, 0, 1),
  c(0, 0, 1, 2, 1, 1)
)
test_loci <- data.frame(
  name = paste0("rs", 1:6),
  chromosome = paste0("chr", c(1, 1, 1, 1, 2, 2)),
  position = as.integer(c(3, 5, 65, 343, 23, 456)),
  genetic_dist = as.double(rep(0, 6)),
  allele_ref = c("A", "T", "C", "G", "C", "T"),
  allele_alt = c("T", "C", "A", "C", "G", "A")
)
test_gt <- gen_tibble(
  x = test_genotypes,
  loci = test_loci,
  indiv_meta = test_indiv_meta,
  quiet = TRUE
)

test_that("genotypes_tbl_class maps the genotypes column to the tibble class", {
  expect_identical(genotypes_tbl_class(test_gt$genotypes), "gen_tbl")
  expect_null(genotypes_tbl_class(1:3))
})

test_that("restore_genotypes_tbl restores the (grouped) gen_tibble class", {
  # a plain tibble
  lost <- as_tibble(test_gt)
  expect_false(inherits(lost, "gen_tbl"))
  expect_identical(class(restore_genotypes_tbl(lost)), class(test_gt))
  # a plain data.frame
  expect_identical(
    class(restore_genotypes_tbl(as.data.frame(test_gt))),
    class(test_gt)
  )
  # a grouped tibble
  test_grouped <- test_gt %>% group_by(population)
  lost_grouped <- lost %>% group_by(population)
  expect_identical(
    class(lost_grouped),
    c("grouped_df", "tbl_df", "tbl", "data.frame")
  )
  restored <- restore_genotypes_tbl(lost_grouped)
  expect_identical(class(restored), class(test_grouped))
  expect_identical(group_data(restored), group_data(test_grouped))
  # a rowwise tibble is ungrouped
  expect_identical(class(restore_genotypes_tbl(rowwise(lost))), class(test_gt))
})

test_that("restore_genotypes_tbl errors on invalid tibbles", {
  expect_error(
    restore_genotypes_tbl(tibble(id = 1:3)),
    "'genotypes' column is missing"
  )
  expect_error(
    restore_genotypes_tbl(tibble(id = 1:3, genotypes = 1:3)),
    "the genotypes column is not of class vctrs_bigSNP"
  )
  # a tibble that already has its class has no method, and would otherwise
  # dispatch to the fallback again
  expect_error(
    restore_genotypes_tbl(test_gt, "my_fun"),
    "my_fun\\(\\) has no method for objects of class 'gen_tbl'"
  )
  expect_error(
    restore_genotypes_tbl(group_by(test_gt, population), "my_fun"),
    "my_fun\\(\\) has no method for objects of class 'grouped_gen_tbl'"
  )
})

test_that("generics restore the class of tibbles that lost it", {
  lost <- as_tibble(test_gt)
  test_grouped <- test_gt %>% group_by(population)
  lost_grouped <- lost %>% group_by(population)
  # gen_tbl methods (via the data.frame fallback)
  expect_identical(qc_report_indiv(lost), qc_report_indiv(test_gt))
  expect_identical(
    qc_report_indiv(as.data.frame(test_gt)),
    qc_report_indiv(test_gt)
  )
  expect_identical(
    pairwise_pop_fst(lost_grouped),
    pairwise_pop_fst(test_grouped)
  )
  expect_identical(
    windows_pop_tajimas_d(lost_grouped, window_size = 3, step_size = 3),
    windows_pop_tajimas_d(test_grouped, window_size = 3, step_size = 3)
  )
  # grouped_gen_tbl methods (via the grouped_df fallback)
  expect_identical(
    loci_alt_freq(lost_grouped),
    loci_alt_freq(test_grouped)
  )
  expect_identical(
    loci_hwe(lost_grouped, type = "matrix"),
    loci_hwe(test_grouped, type = "matrix")
  )
  expect_identical(
    qc_report_loci(lost_grouped),
    qc_report_loci(test_grouped)
  )
  # arguments are passed on
  expect_identical(
    loci_alt_freq(lost_grouped, as_counts = TRUE),
    loci_alt_freq(test_grouped, as_counts = TRUE)
  )
})

test_that("generics work on gen_tibbles returned by pick() and reframe()", {
  test_grouped <- test_gt %>% group_by(population)
  # pick() returns a plain tibble without the grouping variables
  hwe_pick <- test_grouped %>%
    reframe(hwe = loci_hwe(pick(everything())))
  expect_identical(
    hwe_pick$hwe,
    unlist(group_map(test_grouped, ~ loci_hwe(.x)), use.names = FALSE)
  )
  qc_pick <- test_grouped %>%
    summarise(qc = list(qc_report_indiv(pick(everything()))))
  expect_identical(
    qc_pick$qc,
    group_map(test_grouped, ~ qc_report_indiv(.x))
  )
  # reframe() returns a plain tibble, which can be grouped again
  reframed <- test_gt %>%
    reframe(id = id, population = population, genotypes = genotypes)
  expect_false(inherits(reframed, "gen_tbl"))
  expect_identical(
    reframed %>% group_by(population) %>% pairwise_pop_fst(),
    pairwise_pop_fst(test_grouped)
  )
})

test_that("the fallback errors for a class without a method", {
  # a genotypes column of another package, whose tibble class has no method
  registerS3method(
    "genotypes_tbl_class", "fake_genotypes",
    function(genotypes, ...) "fake_tbl",
    envir = asNamespace("tidypopgen")
  )
  fake <- tibble(id = 1:3, genotypes = structure(1:3, class = "fake_genotypes"))
  expect_error(
    pairwise_pop_fst(fake),
    "pairwise_pop_fst\\(\\) has no method for objects of class 'fake_tbl'"
  )
  expect_error(
    loci_hwe(group_by(fake, id)),
    "loci_hwe\\(\\) has no method for objects of class 'grouped_fake_tbl'"
  )
})

test_that("tbl_df methods hand any known genotypes column over", {
  expect_error(
    indiv_missingness(tibble(id = 1:3, genotypes = 1:3)),
    "the genotypes column is not of class vctrs_bigSNP"
  )
  expect_identical(show_loci(as_tibble(test_gt)), show_loci(test_gt))
  expect_identical(count_loci(as_tibble(test_gt)), count_loci(test_gt))
})
