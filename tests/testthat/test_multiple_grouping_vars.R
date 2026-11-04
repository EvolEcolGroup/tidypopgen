# Grouping by several variables must give the same results as grouping by a
# single variable that combines them (with labels pasted with "."), and empty
# groups (which only exist if the user asks for them) must not change the
# results of the other groups.

set.seed(123)
n_indiv <- 24
n_loci <- 40
test_genotypes <- matrix(
  sample(c(0, 1, 2, NA), n_indiv * n_loci,
    replace = TRUE,
    prob = c(0.4, 0.3, 0.25, 0.05)
  ),
  nrow = n_indiv
)
test_indiv_meta <- data.frame(
  id = paste0("ind", seq_len(n_indiv)),
  population = rep(c("pop1", "pop2", "pop3"), each = 8),
  sex = rep(c("F", "M"), times = 12)
)
test_loci <- data.frame(
  name = paste0("rs", seq_len(n_loci)),
  chromosome = rep(c("chr1", "chr2"), each = n_loci / 2),
  position = as.integer(rep(seq_len(n_loci / 2) * 1000, 2)),
  genetic_dist = as.double(rep(0, n_loci)),
  allele_ref = rep("A", n_loci),
  allele_alt = rep("T", n_loci)
)
test_gt <- gen_tibble(
  x = test_genotypes,
  loci = test_loci,
  indiv_meta = test_indiv_meta,
  quiet = TRUE
)
test_gt$pop_sex <- paste(test_gt$population, test_gt$sex, sep = ".")
# grouped by two variables, and by a single variable combining them
gt_two <- test_gt %>% dplyr::group_by(population, sex)
gt_one <- test_gt %>% dplyr::group_by(pop_sex)
group_labels <- c("pop1.F", "pop1.M", "pop2.F", "pop2.M", "pop3.F", "pop3.M")

test_that("loci_* functions work with several grouping variables", {
  loci_funs <- list(
    loci_alt_freq = function(x, ...) loci_alt_freq(x, n_cores = 1, ...),
    loci_maf = function(x, ...) loci_maf(x, n_cores = 1, ...),
    loci_missingness = function(x, ...) loci_missingness(x, n_cores = 1, ...),
    loci_hwe = function(x, ...) loci_hwe(x, n_cores = 1, ...),
    loci_pi = function(x, ...) loci_pi(x, n_cores = 1, ...)
  )
  for (fun in loci_funs) {
    # matrix: same values and labels
    expect_equal(fun(gt_two, type = "matrix"), fun(gt_one, type = "matrix"))
    expect_equal(col_names(fun(gt_two, type = "matrix")), group_labels)
    # list: same values
    expect_equal(fun(gt_two, type = "list"), fun(gt_one, type = "list"))
    # tidy: one column per grouping variable, rather than a single group
    tidy_two <- fun(gt_two, type = "tidy")
    tidy_one <- fun(gt_one, type = "tidy")
    expect_named(tidy_two, c("loci", "population", "sex", "value"))
    expect_equal(tidy_two$loci, tidy_one$loci)
    expect_equal(
      paste(tidy_two$population, tidy_two$sex, sep = "."),
      tidy_one$group
    )
    expect_equal(tidy_two$value, tidy_one$value)
  }
})

test_that("pop_* functions work with several grouping variables", {
  expect_equal(pop_het_obs(gt_two), pop_het_obs(gt_one))
  expect_equal(names(pop_het_obs(gt_two)), group_labels)
  expect_equal(
    pop_het_obs(gt_two, by_locus = TRUE, include_global = TRUE),
    pop_het_obs(gt_one, by_locus = TRUE, include_global = TRUE)
  )
  expect_equal(
    pop_het_exp(gt_two, include_global = TRUE),
    pop_het_exp(gt_one, include_global = TRUE)
  )
  expect_equal(
    pop_fis(gt_two, include_global = TRUE),
    pop_fis(gt_one, include_global = TRUE)
  )
  expect_equal(
    pop_fis(gt_two, method = "WG17", include_global = TRUE),
    pop_fis(gt_one, method = "WG17", include_global = TRUE)
  )
  expect_equal(
    names(pop_fis(gt_two, method = "WG17")),
    group_labels
  )
  expect_equal(
    pop_fst(gt_two, include_global = TRUE),
    pop_fst(gt_one, include_global = TRUE)
  )
  expect_equal(names(pop_fst(gt_two)), group_labels)
  expect_equal(
    pop_global_stats(gt_two, by_locus = TRUE),
    pop_global_stats(gt_one, by_locus = TRUE)
  )
  expect_equal(pop_tajimas_d(gt_two), pop_tajimas_d(gt_one))
})

test_that("pairwise_pop_fst works with several grouping variables", {
  for (method in c("Hudson", "Nei87", "WC84")) {
    tidy_two <- pairwise_pop_fst(gt_two, method = method)
    tidy_one <- pairwise_pop_fst(gt_one, method = method)
    expect_named(
      tidy_two,
      c("population_1", "sex_1", "population_2", "sex_2", "value")
    )
    expect_equal(
      paste(tidy_two$population_1, tidy_two$sex_1, sep = "."),
      tidy_one$pop_sex_1
    )
    expect_equal(
      paste(tidy_two$population_2, tidy_two$sex_2, sep = "."),
      tidy_one$pop_sex_2
    )
    expect_equal(tidy_two$value, tidy_one$value)
    # the pairwise matrix is labelled by the combined labels
    expect_equal(
      pairwise_pop_fst(gt_two, method = method, type = "pairwise"),
      pairwise_pop_fst(gt_one, method = method, type = "pairwise")
    )
    expect_equal(
      pairwise_pop_fst(gt_two,
        method = method, by_locus = TRUE,
        by_locus_type = "matrix"
      )$Fst_by_locus,
      pairwise_pop_fst(gt_one,
        method = method, by_locus = TRUE,
        by_locus_type = "matrix"
      )$Fst_by_locus
    )
  }
  # the tidy output with several grouping variables can be plotted
  expect_s3_class(autoplot(pairwise_pop_fst(gt_two)), "ggplot")
  # and gives the same matrix as with a single grouping variable
  expect_equal(
    ggplot2::layer_data(heatmap_pairwise(pairwise_pop_fst(gt_two))),
    ggplot2::layer_data(heatmap_pairwise(pairwise_pop_fst(gt_one)))
  )
})

test_that("windows and jackknife functions work with several grouping vars", {
  expect_equal(
    windows_pairwise_pop_fst(gt_two,
      window_size = 5, step_size = 5,
      size_unit = "snp"
    ),
    windows_pairwise_pop_fst(gt_one,
      window_size = 5, step_size = 5,
      size_unit = "snp"
    )
  )
  expect_equal(
    jackknife_pairwise_pop_fst(gt_two, block_size = 5, size_unit = "snp"),
    jackknife_pairwise_pop_fst(gt_one, block_size = 5, size_unit = "snp")
  )
  for (type in c("matrix", "tidy")) {
    expect_equal(
      nwise_pop_pbs(gt_two, type = type),
      nwise_pop_pbs(gt_one, type = type)
    )
    expect_equal(
      windows_nwise_pop_pbs(gt_two,
        type = type, window_size = 5,
        step_size = 5, size_unit = "snp"
      ),
      windows_nwise_pop_pbs(gt_one,
        type = type, window_size = 5,
        step_size = 5, size_unit = "snp"
      )
    )
    expect_equal(
      jackknife_nwise_pop_pbs(gt_two,
        type = type, block_size = 5,
        size_unit = "snp"
      ),
      jackknife_nwise_pop_pbs(gt_one,
        type = type, block_size = 5,
        size_unit = "snp"
      )
    )
  }
  for (type in c("matrix", "list")) {
    expect_equal(
      windows_pop_tajimas_d(gt_two,
        type = type, window_size = 5,
        step_size = 5, size_unit = "snp"
      ),
      windows_pop_tajimas_d(gt_one,
        type = type, window_size = 5,
        step_size = 5, size_unit = "snp"
      )
    )
    expect_equal(
      jackknife_pop_tajimas_d(gt_two,
        type = type, block_size = 5,
        size_unit = "snp"
      ),
      jackknife_pop_tajimas_d(gt_one,
        type = type, block_size = 5,
        size_unit = "snp"
      )
    )
  }
  # the tidy outputs have one column per grouping variable
  win_two <- windows_pop_tajimas_d(gt_two,
    type = "tidy", window_size = 5,
    step_size = 5, size_unit = "snp"
  )
  win_one <- windows_pop_tajimas_d(gt_one,
    type = "tidy", window_size = 5,
    step_size = 5, size_unit = "snp"
  )
  expect_named(win_two, c("population", "sex", names(win_one)[-1]))
  expect_equal(paste(win_two$population, win_two$sex, sep = "."), win_one$group)
  expect_equal(win_two[-(1:2)], win_one[-1])
  jack_two <- jackknife_pop_tajimas_d(gt_two,
    type = "tidy", block_size = 5,
    size_unit = "snp"
  )
  jack_one <- jackknife_pop_tajimas_d(gt_one,
    type = "tidy", block_size = 5,
    size_unit = "snp"
  )
  expect_named(jack_two, c("population", "sex", names(jack_one)[-1]))
  expect_equal(jack_two[-(1:2)], jack_one[-1])
})

test_that("qc reports and exports work with several grouping variables", {
  expect_equal(qc_report_loci(gt_two), qc_report_loci(gt_one))
  report_two <- qc_report_indiv(gt_two, kings_threshold = "first")
  report_one <- qc_report_indiv(gt_one, kings_threshold = "first")
  # the group of each individual is a factor, with the levels in the order of
  # the groups
  expect_equal(report_two$group, factor(test_gt$pop_sex, levels = group_labels))
  report_one$group <- factor(report_one$group, levels = group_labels)
  expect_equal(report_two, report_one)
  # the family ID of the plink files are the combined labels
  bed_file <- gt_as_plink(gt_two, file = tempfile(), type = "bed")
  fam <- utils::read.table(sub("\\.bed$", ".fam", bed_file))
  expect_equal(fam$V1, test_gt$pop_sex)
  ped_file <- gt_as_plink(gt_two, file = tempfile(), type = "ped")
  ped <- utils::read.table(ped_file)
  expect_equal(ped$V1, test_gt$pop_sex)
  # and so is the group of a q_matrix
  q_mat <- matrix(stats::runif(n_indiv * 2), ncol = 2)
  q_mat <- q_matrix(q_mat / rowSums(q_mat))
  # (one row per individual and ancestry component)
  expect_equal(
    tidy(q_mat, data = gt_two)$group,
    rep(test_gt$pop_sex, each = 2)
  )
})

# a group emptied with filter(.preserve = TRUE), and the same data without it
gt_empty <- test_gt %>%
  dplyr::group_by(population, sex) %>%
  dplyr::filter(!(population == "pop2" & sex == "F"), .preserve = TRUE)
gt_dropped <- test_gt %>%
  dplyr::group_by(population, sex) %>%
  dplyr::filter(!(population == "pop2" & sex == "F"))
empty_label <- "pop2.F"

test_that("empty groups get NA without changing other groups", {
  expect_equal(dplyr::n_groups(gt_empty), 6)
  expect_equal(dplyr::n_groups(gt_dropped), 5)
  # loci_* functions: the empty group is a column of NA (or NaN)
  for (fun in list(loci_alt_freq, loci_maf, loci_pi)) {
    res <- fun(gt_empty, type = "matrix", n_cores = 1)
    expect_equal(col_names(res), group_labels)
    expect_true(all(is.na(res[, empty_label])))
    expect_equal(
      unname(as.matrix(res[, col_names(res) != empty_label])),
      unname(as.matrix(fun(gt_dropped, type = "matrix", n_cores = 1)))
    )
  }
  # pop_* functions: the empty group is NaN
  for (fun in list(
    function(x) pop_het_obs(x),
    function(x) pop_het_exp(x),
    function(x) pop_fis(x),
    function(x) pop_fis(x, method = "WG17", include_global = TRUE),
    function(x) pop_fst(x, include_global = TRUE)
  )) {
    res <- fun(gt_empty)
    expect_true(is.nan(res[[empty_label]]))
    expect_equal(res[names(res) != empty_label], fun(gt_dropped))
  }
  # pairwise Fst: pairs including the empty group are NaN
  fst_empty <- pairwise_pop_fst(gt_empty)
  fst_dropped <- pairwise_pop_fst(gt_dropped)
  with_empty <- (fst_empty$population_1 == "pop2" & fst_empty$sex_1 == "F") |
    (fst_empty$population_2 == "pop2" & fst_empty$sex_2 == "F")
  expect_true(all(is.nan(fst_empty$value[with_empty])))
  expect_equal(fst_empty$value[!with_empty], fst_dropped$value)
  # qc reports are not affected by the empty group
  expect_equal(qc_report_loci(gt_empty), qc_report_loci(gt_dropped))
  # (but the empty group is still a level of the factor of groups)
  report_empty <- qc_report_indiv(gt_empty, kings_threshold = "first")
  expect_equal(levels(report_empty$group), group_labels)
  report_empty$group <- droplevels(report_empty$group)
  expect_equal(
    report_empty,
    qc_report_indiv(gt_dropped, kings_threshold = "first")
  )
})

test_that("empty groups from .drop = FALSE are kept", {
  # a factor with an unused level
  gt_factor <- test_gt
  gt_factor$pop_factor <- factor(gt_factor$population,
    levels = c("pop1", "pop2", "pop3", "pop4")
  )
  gt_factor <- gt_factor %>% dplyr::group_by(pop_factor, .drop = FALSE)
  gt_pop <- test_gt %>% dplyr::group_by(population)
  res <- loci_missingness(gt_factor, type = "matrix")
  expect_equal(col_names(res), c("pop1", "pop2", "pop3", "pop4"))
  expect_true(all(is.na(res[, "pop4"])))
  expect_equal(
    unname(as.matrix(res[, 1:3])),
    unname(as.matrix(loci_missingness(gt_pop, type = "matrix")))
  )
  expect_equal(
    pop_fst(gt_factor)[1:3],
    pop_fst(gt_pop)
  )
  expect_true(is.nan(pop_fst(gt_factor)[["pop4"]]))
})

test_that("gt_snmf stores the combined group labels", {
  skip_if_not_installed("LEA")
  snmf_res <- gt_snmf(
    gt_two,
    project = "force",
    k = 2,
    n_runs = 1,
    seed = 1
  )
  expect_equal(snmf_res$group, factor(test_gt$pop_sex, levels = group_labels))
})
