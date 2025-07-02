# skip if admixture is not installed
skip_if(
  (system2("which", args = "admixture", stdout = NULL) != 0) &&
    !requireNamespace("tidygenclust", quietly = TRUE)
) # nolint
# set the input file
vcf_path <-
  system.file(
    "/extdata/anolis/punctatus_t70_s10_n46_filtered.recode.vcf.gz",
    package = "tidypopgen"
  )
anole_gt <- gen_tibble(
  vcf_path,
  quiet = TRUE,
  backingfile = tempfile("anolis_")
)
pops_path <- system.file(
  "/extdata/anolis/plot_order_punctatus_n46.csv",
  package = "tidypopgen"
)
pops <- readr::read_csv(pops_path, show_col_types = FALSE)
anole_gt <- anole_gt %>% mutate(id = gsub("punc_", "", .data$id, ))
anole_gt <- anole_gt %>% mutate(population = pops$pop[match(pops$ID, .data$id)])

test_that("gt_admixture error messages", {
  # file that doesn't exist
  invalid_bed <- paste0(tempfile(), ".bed")
  expect_error(
    gt_admixture(
      x = invalid_bed,
      k = 3,
      crossval = FALSE,
      n_cores = 1,
      seed = 123,
      conda_env = "none"
    ),
    "does not exist"
  )
  # random object
  tibble <- tibble::tibble(a = 1:10, b = letters[1:10])
  expect_error(
    gt_admixture(
      x = tibble,
      k = 3,
      crossval = FALSE,
      n_cores = 1,
      seed = 123,
      conda_env = "none"
    ),
    "x must be a gen_tibble or a character"
  )
  # exisitng bed file, but all genotypes missing for one individual
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c"),
    population = c("pop1", "pop1", "pop2")
  )
  test_genotypes <- rbind(
    c(NA, NA, NA, NA, NA, NA),
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
  bed_file <- gt_as_plink(test_gt, file = tempfile(), chromosomes_as_int = TRUE)
  # suppress the warning that admixture has status 255
  # expect printing of error message
  suppressWarnings(expect_error(
    gt_admixture(
      x = bed_file,
      k = 3,
      crossval = FALSE,
      n_cores = 1,
      seed = 123,
      conda_env = "none"
    ),
    "detected that all genotypes are missing for an individual"
  ))
})

test_that("run admixture as single run", {
  # we create a plink file to test the function
  anole_plink <- gt_as_plink(
    anole_gt,
    file = tempfile(),
    chromosomes_as_int = TRUE
  )
  # run admixture
  anole_adm <- gt_admixture(
    anole_plink,
    k = 3,
    crossval = FALSE,
    n_cores = 1,
    seed = 123,
    conda_env = "none"
  )
  # check the output
  expect_true(nrow(anole_adm$Q[[1]]) == nrow(anole_gt))
  expect_true(ncol(anole_adm$Q[[1]]) == 3)
  expect_true(is.null(anole_adm$cv))
  # no create another run and combine them
  anole_adm2 <- gt_admixture(
    anole_plink,
    k = 2,
    crossval = FALSE,
    n_cores = 1,
    seed = 345,
    conda_env = "none"
  )
  anole_adm_comb <- c(anole_adm, anole_adm2)
  expect_true(nrow(anole_adm_comb$Q[[1]]) == nrow(anole_gt))
  expect_true(ncol(anole_adm_comb$Q[[1]]) == 2)
  expect_true(ncol(anole_adm_comb$Q[[2]]) == 3)
  expect_true(all(anole_adm_comb$k == c(2, 3)))
  # run admixture with crossval
  anole_adm3 <- gt_admixture(
    anole_gt,
    k = 3,
    crossval = TRUE,
    n_cores = 2,
    seed = 345,
    conda_env = "none"
  )
  expect_false(is.null(anole_adm3$cv))
  anole_adm4 <- gt_admixture(
    anole_gt,
    k = 3,
    crossval = TRUE,
    n_cores = 2,
    seed = 123,
    conda_env = "none"
  )
  anole_adm_comb2 <- c(anole_adm3, anole_adm4)
  expect_true(ncol(anole_adm_comb2$Q[[1]]) == ncol(anole_adm_comb2$Q[[2]]))
  expect_true(all(anole_adm_comb2$k == c(3, 3)))
  expect_true(all(anole_adm_comb2$cv == c(anole_adm3$cv, anole_adm4$cv)))
})

test_that("run admixture as multiple runs", {
  anole_gt <- anole_gt %>% dplyr::group_by(population)
  anole_adm_cv <- gt_admixture(
    anole_gt,
    k = 2:4,
    n_runs = 2,
    crossval = TRUE,
    n_cores = 2,
    seed = c(123, 234),
    conda_env = "none"
  )
  expect_true(length(anole_adm_cv$k) == 6)
  expect_true(length(anole_adm_cv$Q) == 6)
  expect_true(length(anole_adm_cv$cv) == 6)
  # if we created the object from a grouped gen_tibble,
  # we should have the id and group columns
  expect_true(length(anole_adm_cv$id) == nrow(anole_gt))
  expect_true(length(anole_adm_cv$group) == nrow(anole_gt))
  # error if we have the wrong number of seeds
  expect_error(
    gt_admixture(
      anole_gt,
      k = 2:4,
      n_runs = 2,
      crossval = TRUE,
      n_cores = 2,
      seed = c(123),
      conda_env = "none"
    ),
    "'seed' should be a vector of "
  )
  # plot the crossval
  cross_plot <- autoplot(anole_adm_cv)
  # check that this is indeed a ggplot object
  expect_true(inherits(cross_plot, "ggplot"))
  # plot the admixture results
  expect_error(
    autoplot(anole_adm_cv, type = "barplot"),
    "^You must specify a value for k"
  )
  expect_error(
    autoplot(anole_adm_cv, type = "barplot", k = 2),
    "^You must specify a value for repeat"
  )
  bar_plot <- autoplot(anole_adm_cv, type = "barplot", k = 2, run = 1)
  # check that this is indeed a ggplot object
  expect_true(inherits(bar_plot, "ggplot"))

  # test the reorder of q matrices
  anole_adm_cv_reorder <- gt_admix_reorder_q(anole_adm_cv)
  # check plot ordering
  unord_plot <- autoplot(
    anole_adm_cv,
    type = "barplot",
    k = 3,
    run = 1,
    annotate_group = FALSE
  )
  ord_plot <- autoplot(
    anole_adm_cv,
    type = "barplot",
    k = 3,
    run = 1,
    annotate_group = TRUE
  )
  reord_plot <- autoplot(
    anole_adm_cv_reorder,
    type = "barplot",
    k = 3,
    run = 1,
    annotate_group = TRUE
  )

  expect_identical(ord_plot$data, reord_plot$data)
  expect_false(identical(unord_plot, ord_plot))
  # reordering within plot changes the order further
  ord_within_plot <- autoplot(
    anole_adm_cv_reorder,
    type = "barplot",
    k = 3,
    run = 1,
    annotate_group = TRUE,
    reorder_within_groups = TRUE
  )
  expect_false(identical(ord_plot, ord_within_plot))

  # get q matrix and augment it
  anole_q <- get_q_matrix(anole_adm_cv, k = 3, run = 1)
  augment_q <- augment(anole_q, data = anole_gt)
  expect_true(nrow(augment_q) == nrow(anole_gt))
  # Check that the data are in the same order as the q matrix
  expect_equal(augment_q$group, attr(anole_q, "group"))
  expect_equal(augment_q$id, attr(anole_q, "id"))

  # use the reordered gt_admix object to get the same q matrix
  anole_q_inorder <- get_q_matrix(anole_adm_cv_reorder, k = 3, run = 1)
  # augment with original gen_tibble
  augment_q_inorder <- augment(anole_q_inorder, data = anole_gt)
  # check the resulting augmented gen_tibbles are the same
  expect_equal(augment_q, augment_q_inorder)

  # REINSTATE the test below after fixing tidy and augment

  # tidy matrix with grouped and ungrouped data
  q_tidy_group <- tidy(get_q_matrix(anole_adm_cv, k = 3, run = 1), anole_gt)
  expect_true("group" %in% colnames(q_tidy_group))
  # check that data are in tidy format, each id appears k times
  expect_true(all(table(q_tidy_group$id) == 3))
  # q_tidy_ind <- tidy(get_q_matrix(anole_adm_cv, k=3, run=1), #nolint start
  #                           anole_gt %>% dplyr::ungroup())
  # expect_false("group" %in% colnames(q_tidy_ind)) #nolint end

  # TODO tidy q_matrix with gt that has a different grouping variable
  # to the one used to create the q_matrix and see what happens

  # reorder errors
  # wrong object
  expect_error(
    gt_admix_reorder_q(anole_gt),
    "x must be a gt_admix object"
  )
  # wrong group length
  expect_error(
    gt_admix_reorder_q(anole_adm_cv, group = c(1, 2)),
    paste(
      "The length of the group variable must be the same as the",
      "number of rows in the Q matrix"
    )
  )
  # no group when we have no group information in the gt_admix object
  anole_adm_no_group <- anole_adm_cv
  anole_adm_no_group$group <- NULL
  expect_error(
    gt_admix_reorder_q(anole_adm_no_group),
    paste(
      "You must provide a group variable if there is no grouping",
      "information in the gt_admix object"
    )
  )
  # use gt_adm without id and group
  anole_adm_no_group$id <- NULL
  anole_adm_no_group$group <- NULL
  adm_reord_no_meta <-
    gt_admix_reorder_q(anole_adm_no_group, group = anole_gt$population)
  # now check that we have an id and group elements in the list
  expect_true(length(adm_reord_no_meta$id) == nrow(anole_gt))
  expect_true(length(adm_reord_no_meta$group) == nrow(anole_gt))
})

test_that("assigning factor levels reorders populations in autoplot", {
  # set population to a factor and change the levels
  anole_gt$population <- as.factor(anole_gt$population)
  anole_gt <- anole_gt %>%
    mutate(population = factor(population, levels = c("AF", "Wam", "Eam")))
  anole_gt <- anole_gt %>% group_by(population)
  anole_adm <- gt_admixture(
    anole_gt,
    k = 3,
    crossval = FALSE,
    n_cores = 1,
    seed = 123,
    conda_env = "none"
  )
  expect_error(
    autoplot(anole_adm, type = "cv"),
    "No cross validation error available"
  )
  plt <- autoplot(
    anole_adm,
    type = "barplot",
    k = 3,
    run = 1,
    annotate_group = TRUE,
    arrange_by_group = TRUE,
    arrange_by_indiv = FALSE,
    reorder_within_groups = FALSE
  )
  # and the population order in the plot is changed
  expect_true(all(levels(plt$data$group) == c("AF", "Wam", "Eam")))
  # group of the first individual is therefore "AF"
  expect_true(plt$data$group[1] == "AF")
  # group of the last individual is therefore "Eam"
  expect_true(plt$data$group[138] == "Eam")

  anole_gt <- anole_gt %>% ungroup()
  # reset the population levels
  anole_gt <- anole_gt %>%
    mutate(population = factor(population, levels = c("AF", "Eam", "Wam")))
  anole_gt <- anole_gt %>% group_by(population)

  # reorder the gt_admix object
  anole_adm_reordered1 <-
    gt_admix_reorder_q(anole_adm, group = anole_gt$population)
  plt_b <- autoplot(
    anole_adm_reordered1,
    type = "barplot",
    k = 3,
    run = 1,
    annotate_group = TRUE,
    arrange_by_group = TRUE,
    arrange_by_indiv = FALSE,
    reorder_within_groups = FALSE
  )
  # and the population order changes
  expect_true(all(levels(plt_b$data$group) == c("AF", "Eam", "Wam")))
  # group of the first individual is therefore "AF"
  expect_true(plt_b$data$group[1] == "AF")
  # group of the last individual is therefore "Wam"
  expect_true(plt_b$data$group[138] == "Wam")

  # now test this against the same grouping before running gt_admixture
  anole_gt <- anole_gt %>% ungroup()
  # reset the population levels
  anole_gt <- anole_gt %>%
    mutate(population = factor(population, levels = c("AF", "Eam", "Wam")))
  anole_gt <- anole_gt %>% group_by(population)
  anole_adm2 <- gt_admixture(
    anole_gt,
    k = 3,
    crossval = FALSE,
    n_cores = 1,
    seed = 123,
    conda_env = "none"
  )
  plt2 <- autoplot(
    anole_adm2,
    type = "barplot",
    k = 3,
    run = 1,
    annotate_group = TRUE,
    arrange_by_group = TRUE,
    arrange_by_indiv = FALSE,
    reorder_within_groups = FALSE
  )
  # and the population order changes again
  expect_true(all(levels(plt2$data$group) == c("AF", "Eam", "Wam")))
  # group of the first individual is therefore "AF"
  expect_true(plt2$data$group[1] == "AF")
  # group of the last individual is therefore "Wam"
  expect_true(plt2$data$group[138] == "Wam")

  # both plots should be the same
  expect_identical(plt_b$data, plt2$data)
})

test_that("checking autoplot arrange_ and annotate_ arguments work", {
  # randomly reorder the gt
  order <- sample(seq_len(nrow(anole_gt)))
  anole_gt <- anole_gt %>% arrange(order)
  anole_gt$population <- as.factor(anole_gt$population)
  anole_gt <- anole_gt %>% group_by(population)
  # generate gt_admix object
  anole_adm <- gt_admixture(
    anole_gt,
    k = 3,
    crossval = FALSE,
    n_cores = 1,
    seed = 123,
    conda_env = "none"
  )

  # plot without reordering
  plt1 <- autoplot(
    anole_adm,
    type = "barplot",
    k = 3,
    run = 1,
    annotate_group = FALSE,
    arrange_by_group = FALSE
  )
  expect_true("group" %in% names(plt1$data))
  # check plot data is in the same order as the original gt object
  expect_equal(plt1$data$group, as.factor(rep(anole_gt$population, each = 3)))

  # plot arranged by group
  plt2 <- autoplot(
    anole_adm,
    type = "barplot",
    k = 3,
    run = 1,
    annotate_group = TRUE,
    arrange_by_group = TRUE
  )
  # check the plot has the group variable, arranged with the same levels
  expect_true("group" %in% names(plt2$data))
  expect_equal(levels(anole_gt$population), levels(plt2$data$group))
  # plot arranged by group and individual
  plt3 <- autoplot(
    anole_adm,
    type = "barplot",
    k = 3,
    run = 1,
    annotate_group = TRUE,
    arrange_by_group = TRUE,
    arrange_by_indiv = TRUE,
    reorder_within_groups = TRUE
  )
  # check that the plot has the group variable, arranged with the same levels
  expect_true("group" %in% names(plt3$data))
  expect_equal(levels(anole_gt$population), levels(plt3$data$group))
  # check that the individuals have changed order
  expect_false(identical(plt2$data$id, plt3$data$id))
})


test_that("grouping before running gt_admxiture vs reordering after running gt_admixture", { # nolint
  # Plots grouping BEFORE running admixture
  anole_gt <- anole_gt %>% group_by(population)
  anole_adm_cv <- gt_admixture(
    anole_gt,
    k = 2:4,
    n_runs = 2,
    crossval = TRUE,
    n_cores = 2,
    seed = c(123, 234),
    conda_env = "none"
  )
  basic_plot1 <- autoplot(
    anole_adm_cv,
    type = "barplot",
    k = 3,
    run = 1,
    annotate_group = FALSE,
    arrange_by_group = FALSE,
    arrange_by_indiv = FALSE,
    reorder_within_groups = FALSE
  )
  labels <- autoplot(
    anole_adm_cv,
    type = "barplot",
    k = 3,
    run = 1,
    annotate_group = TRUE,
    arrange_by_group = TRUE,
    arrange_by_indiv = FALSE,
    reorder_within_groups = FALSE
  )

  anole_gt <- anole_gt %>% ungroup()

  # Plots grouping and reordering AFTER running admixture
  anole_adm_cv <- gt_admixture(
    anole_gt,
    k = 2:4,
    n_runs = 2,
    crossval = TRUE,
    n_cores = 2,
    seed = c(123, 234),
    conda_env = "none"
  )
  basic_plot2 <- autoplot(
    anole_adm_cv,
    type = "barplot",
    k = 3,
    run = 1,
    annotate_group = FALSE,
    arrange_by_group = FALSE,
    arrange_by_indiv = FALSE,
    reorder_within_groups = FALSE
  )
  # Expect basic plots to be the same, no reordering has taken place
  expect_identical(basic_plot1$data$id, basic_plot2$data$id)
  expect_identical(basic_plot1$data$percentage, basic_plot2$data$percentage)
  # group and reorder gt_admix
  anole_gt <- anole_gt %>% group_by(population)
  anole_adm_cv_reorder <- gt_admix_reorder_q(
    anole_adm_cv,
    group = anole_gt$population
  )
  # autoplot should now be ordered by group
  # even when arrange_ and annotate_ are false
  reord <- autoplot(
    anole_adm_cv_reorder,
    type = "barplot",
    k = 3,
    run = 1,
    annotate_group = FALSE,
    arrange_by_group = FALSE,
    arrange_by_indiv = FALSE,
    reorder_within_groups = FALSE
  )
  # annotate and arranging by group should keep the same order
  reord_labels <- autoplot(
    anole_adm_cv_reorder,
    type = "barplot",
    k = 3,
    run = 1,
    annotate_group = TRUE,
    arrange_by_group = TRUE,
    arrange_by_indiv = FALSE,
    reorder_within_groups = FALSE
  )
  expect_identical(reord$data, reord_labels$data)

  # the plot from gt_admix taken from grouped gen_tibble should be identical to
  # the plot from a reordered gt_admix object
  expect_identical(labels$data, reord_labels$data)
})
