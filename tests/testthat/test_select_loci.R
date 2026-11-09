test_that("select_loci subsets correctly", {
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c"),
    population = c("pop1", "pop1", "pop2")
  )
  test_genotypes <- rbind(
    c(1, 1, 0, 1, 1, 2),
    c(2, 1, 0, NA, 0, NA),
    c(2, 2, 0, 0, 1, NA)
  )
  test_loci <- data.frame(
    name = c(paste0("rs", 1:4), paste0("x", 1:2)),
    chromosome = c(1, 1, 1, 1, 2, 2),
    position = c(3, 5, 65, 343, 23, 456),
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

  # select snps with an rs
  test_gt_sub <- test_gt %>% select_loci(starts_with("rs"))
  expect_true(!any(c("x1", "x2") %in% loci_names(test_gt_sub)))
  # subsetting by id with reordering
  test_gt_sub <- test_gt %>% select_loci(c(3, 1, 5))
  expect_identical(c("rs3", "rs1", "x1"), loci_names(test_gt_sub))
  # get everything
  test_gt_sub <- test_gt %>% select_loci(everything())
  expect_identical(test_gt, test_gt_sub)
  # use 2:4 range expressions
  test_gt_sub <- test_gt %>% select_loci(2:4)
  expect_identical(test_loci$name[2:4], loci_names(test_gt_sub))
})

test_that("select_loci resolves simple selections as tidyselect does", {
  test_gt <- load_example_gt("gen_tbl")
  n_loci <- count_loci(test_gt)
  # the same gen_tibble with integer loci names
  int_gt <- test_gt
  int_loci <- show_loci(int_gt)
  int_loci$name <- seq_len(n_loci) * 10L
  show_loci(int_gt) <- int_loci
  idx <- c(4, 2, 4)
  # selections that do not depend on the type of names
  selections <- rlang::quos(
    all_of(idx), -all_of(idx), !any_of(c(2, 99)), dplyr::all_of(c(0, 3)),
    -ends_with("2"), c(5, 1), -6
  )
  sel_names <- c("rs3", "rs1")
  char_selections <- rlang::quos(all_of(sel_names), any_of(c("rs2", "zz")))
  int_selections <- rlang::quos(all_of(c("30", "10")), any_of(c("20", "020")))
  cases <- list(
    list(gt = test_gt, sel = c(selections, char_selections)),
    list(gt = int_gt, sel = c(selections, int_selections))
  )
  for (this_case in cases) {
    names_vec <- stats::setNames(seq_len(n_loci), loci_names(this_case$gt))
    for (this_sel in this_case$sel) {
      fast_pos <- select_loci_fast(this_sel, this_case$gt)
      # resolved without tidyselect, with the same result
      expect_false(is.null(fast_pos))
      expect_identical(
        fast_pos,
        unname(tidyselect::eval_select(this_sel, names_vec))
      )
    }
  }
  expect_identical(
    loci_names(select_loci(int_gt, all_of(c("30", "10")))),
    c(30L, 10L)
  )
  # unusual selections fall back on tidyselect, which gives the usual errors
  expect_null(select_loci_fast(rlang::quo(all_of(c(1, NA))), test_gt))
  expect_error(select_loci(test_gt, all_of(99)), "past the end")
  expect_error(select_loci(int_gt, all_of("030")), "don't exist")
  # and combinations are still handled by tidyselect
  expect_identical(
    loci_names(select_loci(test_gt, c(all_of(2), starts_with("rs1")))),
    c("rs2", "rs1")
  )
})

test_that("select_loci evaluates function calls once, as tidyselect does", {
  test_gt <- load_example_gt("gen_tbl")
  n_loci <- count_loci(test_gt)
  # the same gen_tibble with integer loci names
  int_gt <- test_gt
  int_loci <- show_loci(int_gt)
  int_loci$name <- seq_len(n_loci) * 10L
  show_loci(int_gt) <- int_loci
  idx <- c(4, 2, 4)
  # calls (also negated), literals, and injected values
  selections <- rlang::quos(
    sort(idx), -rev(idx), base::rev(c(0, 3)), seq(1, n_loci, by = 2),
    identity(NULL), which(seq_len(n_loci) > 4), 5, -0, !!c(3, 1)
  )
  char_selections <- rlang::quos(
    identity(c("rs3", "rs1")), !setdiff(loci_names(test_gt), "rs2"), "rs2",
    rev(factor(c("rs3", "rs1")))
  )
  int_selections <- rlang::quos(identity(c("30", "10")), "20")
  cases <- list(
    list(gt = test_gt, sel = c(selections, char_selections)),
    list(gt = int_gt, sel = c(selections, int_selections))
  )
  for (this_case in cases) {
    names_vec <- stats::setNames(seq_len(n_loci), loci_names(this_case$gt))
    for (this_sel in this_case$sel) {
      fast_pos <- select_loci_fast(
        select_loci_eval_call(this_sel),
        this_case$gt
      )
      # resolved without tidyselect, with the same result
      expect_false(is.null(fast_pos))
      expect_identical(
        fast_pos,
        unname(tidyselect::eval_select(this_sel, names_vec))
      )
    }
  }
  # the call is evaluated only once, also when its value is left to
  # tidyselect (which gives its usual error)
  n_calls <- 0
  count_calls <- function(value) {
    n_calls <<- n_calls + 1
    value
  }
  expect_identical(
    loci_names(select_loci(test_gt, count_calls(c(3, 1)))),
    c("rs3", "rs1")
  )
  expect_identical(n_calls, 1)
  expect_error(
    select_loci(test_gt, count_calls(c(TRUE, FALSE))),
    "numeric or character"
  )
  expect_identical(n_calls, 2)
  expect_error(select_loci(test_gt, identity(99)), "past the end")
  # calls with selection helpers, and calls that fail, are left to tidyselect
  helper_sel <- rlang::quo(rev(everything()))
  expect_identical(select_loci_eval_call(helper_sel), helper_sel)
  expect_identical(
    loci_names(select_loci(test_gt, rev(everything()))),
    rev(loci_names(test_gt))
  )
  expect_error(select_loci(test_gt, sort(not_a_variable)), "not_a_variable")
})

test_that("select_loci resolves a variable as all_of(), as tidyselect does", {
  test_gt <- load_example_gt("gen_tbl")
  n_loci <- count_loci(test_gt)
  names_vec <- stats::setNames(seq_len(n_loci), loci_names(test_gt))
  idx <- c(4, 2, 4)
  sel_names <- c("rs3", "rs1")
  # the deprecation warning of tidyselect is checked in the next test
  withr::local_options(lifecycle_verbosity = "quiet")
  for (this_sel in rlang::quos(idx, -idx, !sel_names)) {
    fast_pos <- select_loci_fast(this_sel, test_gt)
    # resolved without tidyselect, with the same result
    expect_false(is.null(fast_pos))
    expect_identical(
      fast_pos,
      unname(tidyselect::eval_select(this_sel, names_vec))
    )
  }
  expect_identical(
    select_loci(test_gt, idx),
    select_loci(test_gt, all_of(idx))
  )
  # a locus named like the variable is selected instead, by tidyselect
  rs5 <- 1
  expect_null(select_loci_fast(rlang::quo(rs5), test_gt))
  expect_identical(loci_names(select_loci(test_gt, rs5)), "rs5")
  # values that tidyselect rejects are left to it, to give its error
  keep <- rep(TRUE, n_loci)
  expect_null(select_loci_fast(rlang::quo(keep), test_gt))
  expect_error(select_loci(test_gt, keep), "numeric or character")
})

test_that("select_loci warns about a variable as tidyselect does", {
  test_gt <- load_example_gt("gen_tbl")
  idx <- c(4, 2)
  # show lifecycle warnings every time (not once every 8 hours)
  withr::local_options(lifecycle_verbosity = "warning")
  fast_warning <- expect_warning(
    select_loci(test_gt, idx),
    class = "lifecycle_warning_deprecated"
  )
  tidy_warning <- expect_warning(
    select_loci(test_gt, idx, .engine = "tidyselect"),
    class = "lifecycle_warning_deprecated"
  )
  expect_match(conditionMessage(fast_warning), "external vector")
  expect_identical(
    conditionMessage(fast_warning),
    conditionMessage(tidy_warning)
  )
  # no warning when a locus named like the variable is selected
  rs5 <- 1
  expect_no_warning(select_loci(test_gt, rs5))
})

test_that("select_loci can always use tidyselect", {
  test_gt <- load_example_gt("gen_tbl")
  idx <- c(4, 2, 4)
  # make the fast path fail, to check that it is not used
  local_mocked_bindings(
    select_loci_eval_call = function(...) stop("the fast path was used"),
    select_loci_fast = function(...) stop("the fast path was used")
  )
  expect_identical(
    loci_names(select_loci(test_gt, all_of(idx), .engine = "tidyselect")),
    c("rs4", "rs2")
  )
  expect_identical(
    loci_names(select_loci(test_gt, -sort(idx), .engine = "tidyselect")),
    c("rs1", "rs3", "rs5", "rs6")
  )
  # the fast path is used by default
  expect_error(select_loci(test_gt, all_of(idx)), "the fast path was used")
  expect_error(
    select_loci(test_gt, all_of(idx), .engine = "fast"),
    "should be one of"
  )
})

test_that("select_loci leaves unusual selections to tidyselect", {
  test_gt <- load_example_gt("gen_tbl")
  # no selection, or an empty one, selects no loci
  expect_identical(count_loci(select_loci(test_gt)), 0L)
  expect_identical(count_loci(select_loci(test_gt, c())), 0L)
  # positive and negative positions are combined by tidyselect
  expect_identical(loci_names(select_loci(test_gt, c(1, -2))), "rs1")
  # a variable whose name starts with a dot (like the pronouns of tidyselect)
  .idx <- c(1, 3)
  withr::local_options(lifecycle_verbosity = "quiet")
  expect_identical(loci_names(select_loci(test_gt, .idx)), c("rs1", "rs3"))
  # invalid selections give the errors of tidyselect
  expect_error(select_loci(test_gt, all_of(no_such_var)), "not found")
  expect_error(select_loci(test_gt, starts_with(no_such_var)), "not found")
  expect_error(select_loci(test_gt, no_such_locus), "doesn't exist")
  expect_error(
    select_loci(test_gt, all_of(c("rs1", NA))),
    "missing values"
  )
  expect_error(select_loci(test_gt, all_of(1.5)), "loss of precision")
  expect_error(select_loci(test_gt, all_of(3e9)), "loss of precision")
  expect_error(select_loci(test_gt, c(1, 2.5)), "loss of precision")
})
