test_that("manually resorted table is respected by loci functions", {
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c"),
    population = c("pop1", "pop1", "pop2")
  )
  test_genotypes <- rbind(
    c(1, 1, 0, 1, 1, 0),
    c(2, 1, 0, 0, 0, 0),
    c(2, 2, 0, 0, 1, 1)
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

  # test that loci functions can handle rearranged loci table
  original_freq_alt <- test_gt %>% loci_alt_freq()
  # create new tibble to reorder
  reorder_test_gt <- test_gt
  set.seed(123)
  new_order <- sample(seq_len(nrow(show_loci(test_gt))))
  show_loci(reorder_test_gt) <- show_loci(test_gt)[new_order, ]
  # estimate new frequencies
  reordered_freq_alt <- reorder_test_gt %>% loci_alt_freq()
  # check that the frequencies are the same
  expect_true(all(original_freq_alt[new_order] == reordered_freq_alt))
  # but clumping does generate an error
  expect_no_error(test_gt %>% loci_ld_clump(thr_r2 = 0.2))
  expect_error(
    reorder_test_gt %>% loci_ld_clump(thr_r2 = 0.2),
    "Your loci have been resorted"
  )
})

test_that("gt_update_backingfile correctly updates", {
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c", "d", "e"),
    population = c("pop1", "pop1", "pop2", "pop2", "pop2")
  )
  test_genotypes <- rbind(
    c(1, 1, 0, 1, 1, 0),
    c(2, 1, 0, 0, 0, 0),
    c(2, 2, 0, 0, 1, 1),
    c(NA, 2, 0, NA, NA, 0),
    c(0, 1, 0, 0, 0, 1)
  )
  test_loci <- data.frame(
    name = paste0("rs", 1:6),
    chromosome = paste0("chr", c(2, 1, 1, 1, 1, 2)),
    position = as.integer(c(23, 3, 5, 65, 343, 46)),
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
  subset_reorder_test_gt <- test_gt %>% select_loci(c(2, 4, 5, 1, 6))
  expect_true(is_loci_table_ordered(
    subset_reorder_test_gt,
    error_on_false = FALSE
  ))
  subset_reorder_test_gt <- subset_reorder_test_gt[c(2, 1, 4, 5), ]
  # now save the updated backing matrix
  new_gt <- gt_update_backingfile(subset_reorder_test_gt, quiet = TRUE)

  # check a new .gt is saved
  expect_true(file.exists(bigstatsr::sub_bk(
    gt_get_file_names(new_gt)[2],
    ".gt"
  )))

  # the new gt should be identical to the original one, minus the big indices
  expect_identical(
    show_genotypes(new_gt),
    show_genotypes(subset_reorder_test_gt)
  )
  expect_identical(
    show_loci(new_gt)[, -1],
    show_loci(subset_reorder_test_gt)[, -1]
  )
  # loci big index should now be sequential
  expect_identical(show_loci(new_gt)$big_index, 1:count_loci(new_gt))

  # and if we reload the .gt it is the same
  path_gt <- bigstatsr::sub_bk(gt_get_file_names(new_gt)[2], ".gt")
  reload_gt <- gt_load(path_gt)
  expect_true(all.equal(new_gt, reload_gt, check.attributes = FALSE))

  # if loci are out of order, gt_update_backingfile proceeds, and does not error
  show_loci(subset_reorder_test_gt) <-
    show_loci(subset_reorder_test_gt)[c(1, 3, 2, 4, 5), ]
  out_of_order <- gt_update_backingfile(
    subset_reorder_test_gt,
    quiet = TRUE,
    rm_unsorted_dist = TRUE
  )
  expect_identical(
    show_loci(subset_reorder_test_gt)[, -1],
    show_loci(out_of_order)[, -1]
  )
})

test_that("gt_order_loci reorders and regenerates backingfiles", {
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c"),
    population = c("pop1", "pop1", "pop2")
  )
  test_genotypes <- rbind(
    c(1, 1, 0, 1, 1, 0),
    c(2, 1, 0, 0, 0, 0),
    c(2, 2, 0, 0, 1, 1)
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

  # create new tibble to reorder
  reorder_test_gt <- test_gt
  set.seed(123)
  new_order <- sample(seq_len(nrow(show_loci(test_gt))))
  show_loci(reorder_test_gt) <- show_loci(test_gt)[new_order, ]

  # write out to plink
  path <- paste0(tempfile(), "_outoforder")
  gt_as_plink(reorder_test_gt, path)
  path_bed <- paste0(path, ".bed")

  # read back in
  reorder_test_gt <- gen_tibble(path_bed, quiet = TRUE)
  expect_true(is.integer(show_loci(reorder_test_gt)$position))

  # loci are out of order
  expect_false(is_loci_table_ordered(reorder_test_gt))
  expect_error(
    is_loci_table_ordered(reorder_test_gt, error_on_false = TRUE),
    "Your loci are not sorted within chromosomes"
  )

  # test gt_order_loci
  reorder_test_gt <- gt_order_loci(
    reorder_test_gt,
    use_current_table = FALSE,
    quiet = TRUE
  )

  # loci are now ordered, and backingfiles updated (v2)
  expect_true(is_loci_table_ordered(reorder_test_gt))
  expect_equal(
    normalizePath(gt_get_file_names(reorder_test_gt)),
    normalizePath(c(paste0(path, "_v2.rds"), paste0(path, "_v2.bk")))
  )

  # now save the gen_tibble
  gt_save(reorder_test_gt, quiet = TRUE)

  path_gt <- paste0(path, "_v2.gt")
  # reload and check the backingfile is updated
  reorder_test_gt_reload <- gt_load(path_gt)
  gt_get_file_names(reorder_test_gt_reload)
})

test_that("gt_order_loci catches unsorted and duplicated positions", {
  # gt_order_loci catches unsorted and duplicated positions when
  # use_current_table = TRUE, and duplicated positions
  # when use_current_table = FALSE
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c"),
    population = c("pop1", "pop1", "pop2")
  )
  test_genotypes <- rbind(
    c(1, 1, 0, 1, 1, 0),
    c(2, 1, 0, 0, 0, 0),
    c(2, 2, 0, 0, 1, 1)
  )
  test_loci <- data.frame(
    name = paste0("rs", 1:6),
    chromosome = paste0("chr", c(2, 1, 1, 1, 1, 2)),
    position = as.integer(c(23, 3, 5, 65, 343, 46)),
    genetic_dist = as.double(rep(0, 6)),
    allele_ref = c("A", "T", "C", "G", "C", "T"),
    allele_alt = c("T", "C", NA, "C", "G", "A")
  )
  path <- tempfile()
  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE,
    backingfile = path
  )

  # manually reorder the loci
  loci <- show_loci(test_gt)
  new_order <- c(2, 3, 4, 5, 1, 6)
  show_loci(test_gt) <- show_loci(test_gt)[new_order, ]
  expect_true(is_loci_table_ordered(test_gt))

  # test gt_order_loci use_current_table = TRUE
  test_gt <- gt_order_loci(test_gt, use_current_table = TRUE, quiet = TRUE)
  # now save the gen_tibble
  gt_save(test_gt, quiet = TRUE)

  # loci are still as they were reordered manually
  expect_equal(
    show_loci(test_gt) %>% select(-big_index),
    loci[new_order, ] %>% select(-big_index)
  )

  # but the backingfile has been updated
  expect_equal(
    normalizePath(gt_get_file_names(test_gt)),
    normalizePath(c(paste0(path, "_v2.rds"), paste0(path, "_v2.bk")))
  )

  # if use_current_table = TRUE, but the loci are not ordered, check for errors
  # manually reorder the loci for chromosomes to be out of order
  new_order <- c(5, 1, 2, 3, 4, 6)
  show_loci(test_gt) <- show_loci(test_gt)[new_order, ]
  expect_error(
    gt_order_loci(test_gt, use_current_table = TRUE, quiet = TRUE),
    "All SNPs in a chromosome should be adjacent in the loci table"
  )
  gt_save(test_gt, quiet = TRUE)
  # manually reorder the loci for positions to be out of order
  new_order <- c(3, 4, 2, 5, 1, 6)
  show_loci(test_gt) <- show_loci(test_gt)[new_order, ]
  expect_false(
    is_loci_table_ordered(test_gt),
    "Your loci are not sorted within chromosomes"
  )
  expect_error(
    gt_order_loci(test_gt, use_current_table = TRUE, quiet = TRUE),
    "Your loci are not sorted within chromosomes"
  )
  # if use_current_table = TRUE, but loci contains duplicates, check for errors
  test_loci <- data.frame(
    name = paste0("rs", 1:6),
    chromosome = paste0("chr", c(1, 1, 1, 1, 1, 1)),
    position = as.integer(c(3, 3, 5, 65, 343, 46)),
    genetic_dist = as.double(rep(0, 6)),
    allele_ref = c("A", "T", "C", "G", "C", "T"),
    allele_alt = c("T", "C", NA, "C", "G", "A")
  )
  path <- tempfile()
  expect_warning(test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE,
    allow_duplicates = TRUE,
    backingfile = path
  ), "You have allowed duplicated loci")
  new_order <- c(1, 2, 3, 6, 4, 5)
  # manually reorder so loci are ordered, but have duplicates
  show_loci(test_gt) <- show_loci(test_gt)[new_order, ]
  expect_error(
    is_loci_table_ordered(test_gt, error_on_false = TRUE),
    "Your loci table contains duplicates"
  )

  expect_setequal(
    find_duplicated_loci(test_gt$genotypes,
      error_on_false = FALSE,
      list_duplicates = TRUE
    ),
    c("rs1", "rs2")
  )

  expect_error(
    gt_order_loci(test_gt, use_current_table = TRUE, quiet = TRUE),
    "Your loci table contains duplicates"
  )
  # and if use_current_table = FALSE, check for errors too
  expect_error(
    gt_order_loci(test_gt, use_current_table = FALSE, quiet = TRUE),
    "Your loci table contains duplicates"
  )
})

test_that("gt_order_loci catches unsorted and duplicated genetic_dist", {
  # gt_order_loci catches unsorted and duplicated genetic_dist
  # when use_current_table = TRUE
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c"),
    population = c("pop1", "pop1", "pop2")
  )
  test_genotypes <- rbind(
    c(1, 1, 0, 1, 1, 0),
    c(2, 1, 0, 0, 0, 0),
    c(2, 2, 0, 0, 1, 1)
  )
  test_loci <- data.frame(
    name = paste0("rs", 1:6),
    chromosome = paste0("chr", c(1, 1, 1, 1, 2, 2)),
    position = as.integer(c(3, 5, 65, 343, 27, 83)),
    genetic_dist = as.double(c(0.3, 0.7, 0, 0.2, 0.0, 0.3)),
    allele_ref = c("A", "T", "C", "G", "C", "T"),
    allele_alt = c("T", "C", NA, "C", "G", "A")
  )
  path <- tempfile()
  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE,
    backingfile = path
  )
  expect_true(is_loci_table_ordered(test_gt, ignore_genetic_dist = TRUE))
  expect_false(is_loci_table_ordered(test_gt, ignore_genetic_dist = FALSE))
  # now order it
  ordered_test_gt <-
    gt_order_loci(test_gt,
      use_current_table = FALSE,
      ignore_genetic_dist = TRUE,
      quiet = TRUE
    )

  expect_error(
    gt_order_loci(
      test_gt,
      use_current_table = TRUE,
      ignore_genetic_dist = FALSE,
      quiet = TRUE
    ),
    "Your genetic distances are not sorted within chromosomes"
  )
  expect_error(
    gt_order_loci(
      test_gt,
      use_current_table = FALSE,
      ignore_genetic_dist = FALSE,
      quiet = TRUE
    ),
    "Your genetic distances are not sorted within chromosomes"
  )

  # now we pass the test as we set genetic distances to zero
  expect_true(is_loci_table_ordered(
    ordered_test_gt,
    ignore_genetic_dist = FALSE,
    error_on_false = FALSE
  ))
  # check that genetic_dist is all zero
  expect_equal(show_loci(ordered_test_gt)$genetic_dist, rep(0, 6))
})

test_that("gt_update_backingfile catches unsorted and duplicated genetic_dist", { # nolint
  # gt_update_backingfile catches unsorted and duplicated genetic_dist
  # when rm_unsorted_dist = TRUE
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c"),
    population = c("pop1", "pop1", "pop2")
  )
  test_genotypes <- rbind(
    c(1, 1, 0, 1, 1, 0),
    c(2, 1, 0, 0, 0, 0),
    c(2, 2, 0, 0, 1, 1)
  )
  # Test dist out of order
  test_loci <- data.frame(
    name = paste0("rs", 1:6),
    chromosome = as.character(c(1, 1, 1, 1, 1, 1)),
    position = as.integer(c(3, 5, 25, 46, 65, 343)),
    genetic_dist = c(0.1, 0.3, 0.2, 0.4, 0.5, 0.6),
    allele_ref = c("A", "T", "C", "G", "C", "T"),
    allele_alt = c("T", "C", NA, "C", "G", "A")
  )
  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE,
    backingfile = tempfile()
  )
  expect_equal(show_loci(test_gt)$genetic_dist, c(0.1, 0.3, 0.2, 0.4, 0.5, 0.6))

  test_gt <- gt_update_backingfile(test_gt,
    rm_unsorted_dist = TRUE,
    quiet = TRUE
  )
  # after updating backingfile, genetic_dist should all be zero
  expect_equal(show_loci(test_gt)$genetic_dist, c(0, 0, 0, 0, 0, 0))
})


test_that("is_loci_table_ordered catches unsorted genetic_dist", { # nolint
  # is_loci_table_ordered catches unsorted genetic_dist
  # when ignore_genetic_dist = FALSE
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c"),
    population = c("pop1", "pop1", "pop2")
  )
  test_genotypes <- rbind(
    c(1, 1, 0, 1, 1, 0),
    c(2, 1, 0, 0, 0, 0),
    c(2, 2, 0, 0, 1, 1)
  )
  test_loci <- data.frame(
    name = paste0("rs", 1:6),
    chromosome = as.character(c(1, 1, 1, 1, 1, 1)),
    position = as.integer(c(3, 5, 25, 46, 65, 343)),
    genetic_dist = c(0.1, 0.3, 0.2, 0.4, 0.5, 0.6),
    allele_ref = c("A", "T", "C", "G", "C", "T"),
    allele_alt = c("T", "C", NA, "C", "G", "A")
  )
  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE,
    backingfile = tempfile()
  )
  expect_false(is_loci_table_ordered(test_gt, ignore_genetic_dist = FALSE))
  expect_error(
    is_loci_table_ordered(
      test_gt,
      error_on_false = TRUE,
      ignore_genetic_dist = FALSE
    ),
    "Your genetic distances are not sorted within chromosomes"
  )
  # is_loci_table_ordered ignores duplication and returns TRUE warning
  # when all genetic_dist is set to 0
  test_loci <- data.frame(
    name = paste0("rs", 1:6),
    chromosome = as.character(c(1, 1, 1, 1, 1, 1)),
    position = as.integer(c(3, 4, 25, 46, 65, 343)),
    genetic_dist = as.double(rep(0, 6)),
    allele_ref = c("A", "T", "C", "G", "C", "T"),
    allele_alt = c("T", "C", NA, "C", "G", "A")
  )
  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE,
    backingfile = tempfile()
  )
  expect_true(is_loci_table_ordered(test_gt, ignore_genetic_dist = FALSE))
})

test_that("is_loci_table_ordered catches unsorted and duplicated positions", {
  # is_loci_table_ordered catches duplicated positions
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c"),
    population = c("pop1", "pop1", "pop2")
  )
  test_genotypes <- rbind(
    c(1, 1, 0, 1, 1, 0),
    c(2, 1, 0, 0, 0, 0),
    c(2, 2, 0, 0, 1, 1)
  )
  test_loci <- data.frame(
    name = paste0("rs", 1:6),
    chromosome = as.character(c(1, 1, 1, 1, 1, 1)),
    position = as.integer(c(3, 3, 25, 46, 65, 343)),
    genetic_dist = as.double(rep(0, 6)),
    allele_ref = c("A", "T", "C", "G", "C", "T"),
    allele_alt = c("T", "C", NA, "C", "G", "A")
  )
  expect_warning(test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE,
    allow_duplicates = TRUE,
    backingfile = tempfile()
  ), "You have allowed duplicated loci")
  expect_false(is_loci_table_ordered(test_gt))
  expect_error(
    is_loci_table_ordered(test_gt, error_on_false = TRUE),
    "Your loci table contains duplicates"
  )
  test_loci <- data.frame(
    name = paste0("rs", 1:6),
    chromosome = paste0("chr", c(2, 1, 1, 1, 1, 2)),
    position = as.integer(c(23, 3, 5, 65, 343, 46)),
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
  # subset and reorder
  expect_false(is_loci_table_ordered(test_gt, error_on_false = FALSE))
  expect_error(
    is_loci_table_ordered(test_gt, error_on_false = TRUE),
    "All SNPs in a chromosome should be adjacent"
  )
  test_loci <- data.frame(
    name = paste0("rs", 1:6),
    chromosome = paste0("chr", c(1, 1, 1, 1, 2, 2)),
    position = as.integer(c(5, 3, 65, 343, 23, 456)),
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
  # loci are out of order
  expect_false(is_loci_table_ordered(test_gt))
  expect_error(
    is_loci_table_ordered(test_gt, error_on_false = TRUE),
    "Your loci are not sorted within chromosomes"
  )
})


test_that("check updated positions/distances are inherited by the merged gt", {
  raw_path_pop_b <- system.file("extdata/pop_b.bed", package = "tidypopgen")
  bigsnp_path_b <- bigsnpr::snp_readBed(
    raw_path_pop_b,
    backingfile = tempfile("test_b_")
  )
  pop_b_gt <- gen_tibble(bigsnp_path_b, quiet = TRUE)
  raw_path_pop_a <- system.file("extdata/pop_a.bed", package = "tidypopgen")
  bigsnp_path_a <- bigsnpr::snp_readBed(
    raw_path_pop_a,
    backingfile = tempfile("test_a_")
  )
  pop_a_gt <- gen_tibble(bigsnp_path_a, quiet = TRUE)
  # Edit genetic_dist in loci table of pop_a_gt
  show_loci(pop_a_gt)$genetic_dist <- c(seq(0.01, 0.16, 0.01))
  # Update backingfiles to store this change
  pop_a_gt <- gt_update_backingfile(
    pop_a_gt,
    rm_unsorted_dist = FALSE,
    quiet = TRUE
  )
  # check we're now using the correct backingfiles
  # gt_get_file_names(pop_a_gt) #nolint
  merged_gen <- rbind.gen_tbl(
    pop_a_gt,
    pop_b_gt,
    flip_strand = TRUE,
    quiet = TRUE,
    backingfile = tempfile()
  )
  # check that the genetic_dist is as expected in the merged gt
  expect_equal(
    show_loci(merged_gen)$genetic_dist,
    c(seq(0.02, 0.06, 0.01), 0.08, 0.09, seq(0.12, 0.16, 0.01))
  )
})


test_that("cast_chromosome_to_int transforms correctly from character", {
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c"),
    population = c("pop1", "pop1", "pop2")
  )
  test_genotypes <- rbind(
    c(1, 1, 0, 1, 1, 0),
    c(2, 1, 0, 0, 0, 0),
    c(2, 2, 0, 0, 1, 1)
  )
  test_loci <- data.frame(
    name = paste0("rs", 1:6),
    chromosome = as.character(c(2, 13, 8, 1, 2, 3)),
    position = as.integer(c(23, 3, 5, 65, 343, 46)),
    genetic_dist = as.double(rep(0, 6)),
    allele_ref = c("A", "T", "C", "G", "C", "T"),
    allele_alt = c("T", "C", NA, "C", "G", "A")
  )
  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE,
    backingfile = tempfile()
  )
  expect_equal(
    cast_chromosome_to_int(show_loci(test_gt)$chromosome),
    c(2, 13, 8, 1, 2, 3)
  )
  # use gt_order_loci to reorder
  test_gt <- gt_order_loci(test_gt, use_current_table = FALSE, quiet = TRUE)
  expect_equal(
    cast_chromosome_to_int(show_loci(test_gt)$chromosome),
    c(1, 2, 2, 3, 8, 13)
  )
})


test_that("cast_chromosome_to_int transforms correctly from factor", {
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c"),
    population = c("pop1", "pop1", "pop2")
  )
  test_genotypes <- rbind(
    c(1, 1, 0, 1, 1, 0),
    c(2, 1, 0, 0, 0, 0),
    c(2, 2, 0, 0, 1, 1)
  )
  test_loci <- data.frame(
    name = paste0("rs", 1:6),
    chromosome = as.factor(c(2, 13, 8, 1, 2, 3)),
    position = as.integer(c(23, 3, 5, 65, 343, 46)),
    genetic_dist = as.double(rep(0, 6)),
    allele_ref = c("A", "T", "C", "G", "C", "T"),
    allele_alt = c("T", "C", NA, "C", "G", "A")
  )
  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE,
    backingfile = tempfile()
  )
  # check chr_int makes sense
  expect_equal(
    cast_chromosome_to_int(show_loci(test_gt)$chromosome),
    c(2, 13, 8, 1, 2, 3)
  )
  # use gt_order_loci to reorder
  test_gt <- gt_order_loci(test_gt, use_current_table = FALSE, quiet = TRUE)
  expect_equal(
    cast_chromosome_to_int(show_loci(test_gt)$chromosome),
    c(1, 2, 2, 3, 8, 13)
  )
  expect_equal(
    as.character(show_loci(test_gt)$chromosome),
    as.character(c(1, 2, 2, 3, 8, 13))
  )
})
