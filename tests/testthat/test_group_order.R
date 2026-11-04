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

test_that("tidy, augment and autoplot of a q_matrix keep the order of groups", {
  q <- seq(0.1, 0.7, length.out = nrow(test_gt))
  q_mat <- q_matrix(matrix(c(q, 1 - q), ncol = 2))
  factor_gt <- test_gt %>% dplyr::group_by(pop_factor)
  number_gt <- test_gt %>% dplyr::group_by(pop_number)
  # the groups from grouped data keep their type
  expect_equal(augment(q_mat, data = factor_gt)$group, test_gt$pop_factor)
  expect_equal(
    tidy(q_mat, data = number_gt)$group,
    rep(test_gt$pop_number, each = 2)
  )
  # autoplot() orders the groups as the grouped functions do, whether they
  # come from the data or from the attributes of the q_matrix
  plot_labels <- function(plt) ggplot2::layer_scales(plt)$x$labels
  expect_equal(
    plot_labels(autoplot(q_mat, data = factor_gt)),
    c("pop3", "pop1", "pop2")
  )
  expect_equal(
    plot_labels(autoplot(q_mat, data = number_gt)),
    c("1", "2", "10")
  )
  q_attr <- q_mat
  attr(q_attr, "id") <- test_gt$id
  attr(q_attr, "group") <- test_gt$pop_factor
  expect_equal(plot_labels(autoplot(q_attr)), c("pop3", "pop1", "pop2"))
  # a q_matrix with ids but no groups (e.g. from an ungrouped gen_tibble)
  # takes the groups from grouped data, matched by id
  q_id <- q_mat
  attr(q_id, "id") <- test_gt$id
  rev_gt <- test_gt %>%
    dplyr::arrange(dplyr::desc(id)) %>%
    dplyr::group_by(pop_factor)
  expect_equal(
    tidy(q_id, data = rev_gt)$group,
    rep(test_gt$pop_factor, each = 2)
  )
  augment_id <- augment(q_id, data = rev_gt)
  expect_equal(augment_id$group, rev_gt$pop_factor)
  expect_equal(augment_id$.Q1, q[match(rev_gt$id, test_gt$id)])
  expect_equal(
    plot_labels(autoplot(q_id, data = rev_gt)),
    c("pop3", "pop1", "pop2")
  )
})

test_that("reorder_within_groups keeps the order of the groups", {
  q <- seq(0.1, 0.7, length.out = nrow(test_gt))
  q_attr <- q_matrix(matrix(c(q, 1 - q), ncol = 2))
  attr(q_attr, "id") <- test_gt$id
  attr(q_attr, "group") <- test_gt$pop_factor
  plot_labels <- function(plt) ggplot2::layer_scales(plt)$x$labels
  # arranged by group: in the order of the levels
  expect_equal(
    plot_labels(autoplot(q_attr, reorder_within_groups = TRUE)),
    c("pop3", "pop1", "pop2")
  )
  # not arranged by group: in their current order (of first appearance), as
  # without reorder_within_groups
  expect_equal(
    plot_labels(autoplot(q_attr, arrange_by_group = FALSE)),
    c("pop1", "pop2", "pop3")
  )
  expect_equal(
    plot_labels(autoplot(q_attr,
      arrange_by_group = FALSE,
      reorder_within_groups = TRUE
    )),
    c("pop1", "pop2", "pop3")
  )
})

test_that("heatmap_pairwise keeps the order of a tidy input", {
  axis_labels <- function(plt) ggplot2::layer_scales(plt)$x$labels
  # groups, in the order of the levels
  factor_gt <- test_gt %>% dplyr::group_by(pop_factor)
  fst_mat <- pairwise_pop_fst(factor_gt, type = "pairwise", method = "Nei87")
  fst_tidy <- pairwise_pop_fst(factor_gt, type = "tidy", method = "Nei87")
  expect_equal(
    axis_labels(heatmap_pairwise(fst_mat)),
    c("pop3", "pop1", "pop2")
  )
  expect_equal(
    axis_labels(heatmap_pairwise(fst_tidy)),
    c("pop3", "pop1", "pop2")
  )
  # individuals, in the order of the gen_tibble
  rev_gt <- test_gt %>% dplyr::arrange(dplyr::desc(id))
  king_mat <- pairwise_king(rev_gt, as_matrix = TRUE)
  king_tidy <- pairwise_king(rev_gt)
  expect_equal(axis_labels(heatmap_pairwise(king_mat)), rev_gt$id)
  expect_equal(axis_labels(heatmap_pairwise(king_tidy)), rev_gt$id)
})
