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
