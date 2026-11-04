# Groups are in the order of dplyr::group_keys() (the levels of a factor,
# numbers in numeric order), and outputs keep that order when they are
# sorted or plotted.

test_gt <- load_example_gt("gen_tbl")
# a factor with levels that are not in alphabetical order, and an empty level
test_gt$pop_factor <- factor(test_gt$population,
  levels = c("pop3", "pop_empty", "pop1", "pop2")
)
# numbers, which would sort as "10" < "2" as characters
test_gt$pop_number <- c(pop1 = 10, pop2 = 2, pop3 = 1)[test_gt$population]

test_that("tidy outputs keep the type of the grouping variable", {
  factor_gt <- test_gt %>% dplyr::group_by(pop_factor, .drop = FALSE)
  number_gt <- test_gt %>% dplyr::group_by(pop_number)
  tidy_funs <- list(
    function(x) loci_alt_freq(x, type = "tidy"),
    function(x) loci_maf(x, type = "tidy"),
    function(x) loci_missingness(x, type = "tidy"),
    function(x) loci_hwe(x, type = "tidy"),
    function(x) loci_pi(x, type = "tidy"),
    function(x) {
      windows_pop_tajimas_d(x,
        type = "tidy", window_size = 3, step_size = 2,
        size_unit = "snp"
      )
    },
    function(x) {
      jackknife_pop_tajimas_d(x,
        type = "tidy", block_size = 2,
        size_unit = "snp"
      )
    }
  )
  for (fun in tidy_funs) {
    # a factor keeps its levels, and the rows are in their order
    res <- fun(factor_gt)
    expect_s3_class(res$group, "factor")
    expect_equal(levels(res$group), levels(test_gt$pop_factor))
    expect_equal(
      as.character(unique(res$group)),
      levels(test_gt$pop_factor)
    )
    # numbers stay numbers, in numeric order
    res <- fun(number_gt)
    expect_type(res$group, "double")
    expect_equal(unique(res$group), c(1, 2, 10))
  }
  # matrix outputs are labelled with the groups, in the same order
  expect_equal(
    colnames(loci_alt_freq(number_gt, type = "matrix")),
    c("1", "2", "10")
  )
  expect_equal(
    colnames(loci_alt_freq(factor_gt, type = "matrix")),
    levels(test_gt$pop_factor)
  )
})

test_that("groups of individuals from several grouping variables are ordered", {
  two_gt <- test_gt %>% dplyr::group_by(pop_factor, pop_number)
  # the labels, in the order of the groups (alphabetically, pop3.1 would be
  # last)
  labels <- c("pop3.1", "pop1.10", "pop2.2")
  expect_equal(
    paste(dplyr::group_keys(two_gt)$pop_factor,
      dplyr::group_keys(two_gt)$pop_number,
      sep = "."
    ),
    labels
  )
  # qc_report_indiv() gives a factor, with levels in the order of the groups
  report <- qc_report_indiv(two_gt)
  expect_s3_class(report$group, "factor")
  expect_equal(levels(report$group), labels)
  expect_equal(
    as.character(report$group),
    paste(test_gt$pop_factor, test_gt$pop_number, sep = ".")
  )
  # so does the group stored by gt_admixture() and gt_snmf(), and
  # gt_admix_reorder_q() then sorts the individuals in the order of the groups
  q <- seq(0.1, 0.7, length.out = nrow(test_gt))
  admix_obj <- structure(
    list(
      k = 2L,
      Q = list(q_matrix(matrix(c(q, 1 - q), ncol = 2))),
      id = test_gt$id,
      group = tidypopgen:::.group_of_rows(two_gt),
      algorithm = "ADMIXTURE"
    ),
    class = c("gt_admix", "list")
  )
  expect_equal(levels(admix_obj$group), labels)
  reordered <- gt_admix_reorder_q(admix_obj)
  expect_equal(as.character(unique(reordered$group)), labels)
  expect_equal(
    reordered$id,
    test_gt$id[order(match(test_gt$population, c("pop3", "pop1", "pop2")))]
  )
})
