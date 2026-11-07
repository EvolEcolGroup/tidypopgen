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

test_that("the jackknife generics restore the class of tibbles that lost it", {
  lost <- as_tibble(test_gt)
  test_grouped <- test_gt %>% group_by(population)
  lost_grouped <- lost %>% group_by(population)
  # gen_tbl methods (via the data.frame fallback)
  expect_identical(
    jackknife_pop_tajimas_d(lost, block_size = 2),
    jackknife_pop_tajimas_d(test_gt, block_size = 2)
  )
  expect_identical(
    jackknife_pop_tajimas_d(lost_grouped, type = "tidy", block_size = 2),
    jackknife_pop_tajimas_d(test_grouped, type = "tidy", block_size = 2)
  )
  expect_identical(
    jackknife_pairwise_pop_fst(lost_grouped, block_size = 2),
    jackknife_pairwise_pop_fst(test_grouped, block_size = 2)
  )
  # PBS needs three populations
  test_trio <- test_gt %>% mutate(trio = rep(c("a", "b", "c"), each = 2))
  expect_identical(
    jackknife_nwise_pop_pbs(
      as_tibble(test_trio) %>% group_by(trio),
      block_size = 2
    ),
    jackknife_nwise_pop_pbs(test_trio %>% group_by(trio), block_size = 2)
  )
  # the arguments are checked by the methods
  expect_error(
    jackknife_pop_tajimas_d(lost_grouped, block_size = 2, n_prj = 4),
    "must be empty"
  )
  expect_error(
    jackknife_pairwise_pop_fst(test_grouped, block_size = 2, n_prj = 4),
    "must be empty"
  )
})

test_that("the SFS generics restore the class of tibbles that lost it", {
  lost <- as_tibble(test_gt)
  test_grouped <- test_gt %>% group_by(population)
  lost_grouped <- lost %>% group_by(population)
  # tbl_df methods hand the genotypes column over
  expect_identical(pop_sfs(lost, n_proj = 4), pop_sfs(test_gt, n_proj = 4))
  expect_identical(pop_sfs_preview(lost), pop_sfs_preview(test_gt))
  # grouped_gen_tbl methods (via the grouped_df fallback)
  expect_identical(
    pop_sfs(lost_grouped, n_proj = 4),
    pop_sfs(test_grouped, n_proj = 4)
  )
  expect_identical(
    pop_sfs_preview(lost_grouped),
    pop_sfs_preview(test_grouped)
  )
  # gen_tbl methods (via the data.frame fallback)
  expect_identical(
    pairwise_pop_sfs(lost_grouped, n_proj = 4),
    pairwise_pop_sfs(test_grouped, n_proj = 4)
  )
  expect_identical(
    nwise_pop_sfs(lost_grouped, n_proj = 4, type = "array"),
    nwise_pop_sfs(test_grouped, n_proj = 4, type = "array")
  )
  lost_dir <- withr::local_tempdir()
  grouped_dir <- withr::local_tempdir()
  write_sfs(lost_grouped, out_dir = lost_dir, n_proj = 4)
  write_sfs(test_grouped, out_dir = grouped_dir, n_proj = 4)
  files <- list.files(grouped_dir, recursive = TRUE)
  expect_setequal(files, c("pop1-4.sfs", "pop2-4.sfs", "pop1-pop2.sfs"))
  expect_identical(list.files(lost_dir, recursive = TRUE), files)
  for (f in files) {
    expect_identical(
      readLines(file.path(lost_dir, f)),
      readLines(file.path(grouped_dir, f))
    )
  }
  # a data frame without genotypes (e.g. a tidy spectrum) is not taken for a
  # gen_tibble that lost its class
  expect_error(
    write_sfs(pop_sfs(test_grouped, n_proj = 4), tempfile()),
    "should be an sfs object"
  )
  # the arguments are checked by the methods
  expect_error(pop_sfs(lost_grouped, n_prj = 4), "must be empty")
  expect_error(pairwise_pop_sfs(lost_grouped, n_prj = 4), "must be empty")
  expect_error(nwise_pop_sfs(test_grouped, n_prj = 4), "must be empty")
  # bootstrap functions (gen_tbl methods, via the data.frame fallback)
  expect_identical(
    withr::with_seed(1, bootstrap_pop_sfs(lost_grouped,
      n_bootstrap = 2, block_size = 2, n_proj = 4
    )),
    withr::with_seed(1, bootstrap_pop_sfs(test_grouped,
      n_bootstrap = 2, block_size = 2, n_proj = 4
    ))
  )
  expect_identical(
    withr::with_seed(1, bootstrap_pop_sfs(lost,
      n_bootstrap = 2, block_size = 2, n_proj = 4
    )),
    withr::with_seed(1, bootstrap_pop_sfs(test_gt,
      n_bootstrap = 2, block_size = 2, n_proj = 4
    ))
  )
  expect_identical(
    withr::with_seed(1, bootstrap_pairwise_pop_sfs(lost_grouped,
      n_bootstrap = 2, block_size = 2, n_proj = 4
    )),
    withr::with_seed(1, bootstrap_pairwise_pop_sfs(test_grouped,
      n_bootstrap = 2, block_size = 2, n_proj = 4
    ))
  )
  expect_identical(
    withr::with_seed(1, bootstrap_nwise_pop_sfs(lost_grouped,
      n_bootstrap = 2, block_size = 2, n_proj = 4, type = "array"
    )),
    withr::with_seed(1, bootstrap_nwise_pop_sfs(test_grouped,
      n_bootstrap = 2, block_size = 2, n_proj = 4, type = "array"
    ))
  )
  expect_error(
    bootstrap_pop_sfs(lost_grouped, block_size = 2, n_prj = 4),
    "must be empty"
  )
  # pick() returns a plain tibble without the grouping variables
  sfs_pick <- test_grouped %>%
    summarise(sfs = list(
      pop_sfs(pick(everything()), n_proj = 4, type = "array")
    ))
  expect_identical(
    sfs_pick$sfs,
    group_map(test_grouped, ~ pop_sfs(.x, n_proj = 4, type = "array"))
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
  expect_error(
    pairwise_pop_sfs(fake),
    "pairwise_pop_sfs\\(\\) has no method for objects of class 'fake_tbl'"
  )
  expect_error(
    pop_sfs(group_by(fake, id)),
    "pop_sfs\\(\\) has no method for objects of class 'grouped_fake_tbl'"
  )
})

test_that("tbl_df methods hand any known genotypes column over", {
  expect_error(
    indiv_missingness(tibble(id = 1:3, genotypes = 1:3)),
    "the genotypes column is not of class vctrs_bigSNP"
  )
  expect_error(
    pop_sfs(tibble(id = 1:3, genotypes = 1:3)),
    "the genotypes column is not of class vctrs_bigSNP"
  )
  expect_identical(show_loci(as_tibble(test_gt)), show_loci(test_gt))
  expect_identical(count_loci(as_tibble(test_gt)), count_loci(test_gt))
})
