
# skip if admixture is not installed
skip_if((system2("which", args = "admixture", stdout = NULL) != 0) && !requireNamespace("fastmixturer", quietly = TRUE))
# set the input file
vcf_path <- system.file("/extdata/anolis/punctatus_t70_s10_n46_filtered.recode.vcf.gz",
                        package = "tidypopgen")
anole_gt <- gen_tibble(vcf_path, quiet = TRUE, backingfile = tempfile("anolis_"))
pops_path <- system.file("/extdata/anolis/plot_order_punctatus_n46.csv",
                         package = "tidypopgen")
pops <- readr::read_csv(pops_path, show_col_types = FALSE)
anole_gt <- anole_gt %>% mutate(id = gsub('punc_',"",.data$id,))
anole_gt <- anole_gt %>% mutate(population = pops$pop[match(pops$ID,.data$id)])

test_that("run admixture as single run", {
  # we create a plink file to test the function
  anole_plink <- gt_as_plink(anole_gt, file = tempfile(), chromosomes_as_int=TRUE)
  # run admixture
  anole_adm <- gt_admixture(anole_plink, k = 3, crossval = FALSE, n_cores = 1, seed = 123, conda_env = "none")
  # check the output
  expect_true(nrow(anole_adm$Q[[1]])==nrow(anole_gt))
  expect_true(ncol(anole_adm$Q[[1]])==3)
  expect_true(is.null(anole_adm$cv))
  # no create another run and combine them
  anole_adm2 <- gt_admixture(anole_plink, k = 2, crossval = FALSE, n_cores = 1, seed = 345, conda_env = "none")
  anole_adm_comb <- c(anole_adm, anole_adm2)
  expect_true(nrow(anole_adm_comb$Q[[1]])==nrow(anole_gt))
  expect_true(ncol(anole_adm_comb$Q[[1]])==3)
  expect_true(ncol(anole_adm_comb$Q[[2]])==2)
  expect_true(all(anole_adm_comb$k==c(3,2)))
  # run admixture with crossval
  anole_adm3 <- gt_admixture(anole_gt, k = 3, crossval = TRUE, n_cores = 2, seed = 345, conda_env = "none")
  expect_false(is.null(anole_adm3$cv))
  anole_adm4 <- gt_admixture(anole_gt, k = 3, crossval = TRUE, n_cores = 2, seed = 123, conda_env = "none")
  anole_adm_comb2 <- c(anole_adm3, anole_adm4)
  # TODO write some check for the object above

})

test_that("run admixture as multiple runs", {
  anole_gt <-  anole_gt %>% dplyr::group_by(population)
  anole_adm_cv <- gt_admixture(anole_gt, k = 2:4, n_runs =2, crossval = TRUE, n_cores = 2, seed = c(123,234), conda_env = "none")
  expect_true(length(anole_adm_cv$k)==6)
  expect_true(length(anole_adm_cv$Q)==6)
  expect_true(length(anole_adm_cv$cv)==6)
  # if we created the object from a grouped gen_tibble, we should have the id and group columns
  expect_true(length(anole_adm_cv$id)==nrow(anole_gt))
  expect_true(length(anole_adm_cv$group)==nrow(anole_gt))
  # error if we have the wrong number of seeds
  expect_error(gt_admixture(anole_gt, k = 2:4, n_runs =2, crossval = TRUE, n_cores = 2, seed = c(123), conda_env = "none"),
               "'seeds' should be a vector of ")
  # plot the crossval
  cross_plot <- autoplot(anole_adm_cv)
  # check that this is indeed a ggplot object
  expect_true(inherits(cross_plot, "ggplot"))
  # plot the admixture results
  expect_error(autoplot(anole_adm_cv, type = "barplot"),
               "^You must specify a value for k")
  expect_error(autoplot(anole_adm_cv, type = "barplot", k=2),
               "^You must specify a value for repeat")
  bar_plot <- autoplot(anole_adm_cv, type = "barplot", k=2,run = 1)
  # check that this is indeed a ggplot object
  expect_true(inherits(bar_plot, "ggplot"))

  # test the reorder of q matrices
  anole_adm_cv_reorder <- gt_admix_reorder_q(anole_adm_cv)
  # check plot ordering
  unord_plot <- autoplot(anole_adm_cv, type = "barplot", k=3,run = 1, annotate_group=FALSE)
  ord_plot <- autoplot(anole_adm_cv_reorder, type = "barplot", k=3,run = 1, annotate_group=TRUE)
  reord_plot <- autoplot(anole_adm_cv_reorder, type = "barplot", k=3,run = 1, annotate_group=TRUE)
  expect_identical (ord_plot, reord_plot)
  expect_false (identical (unord_plot, ord_plot))
  # reordering within plot changes the order further
  ord_within_plot <- autoplot(anole_adm_cv_reorder, type = "barplot", k=3,run = 1,
                              annotate_group=TRUE, reorder_within_groups=TRUE)
  expect_false (identical (ord_plot, ord_within_plot))

  # get q matrix and augment it
  anole_q <- get_q_matrix(anole_adm_cv, k=3, run=1)
  augment_q <- augment(anole_q, data= anole_gt)
  expect_true(nrow(augment_q)==nrow(anole_gt))
  # TODO should we try to check that the data are in the same order as the q matrix?

  # reorder errors
  # wrong object
  expect_error(gt_admix_reorder_q(anole_gt),
               "x must be a gt_admix object")
  # wrong group length
  expect_error(gt_admix_reorder_q(anole_adm_cv, group=c(1,2)),
               "The length of the group variable must be the same as the number of rows in the Q matrix")
  # no group when we have no group info
  anole_adm_no_group <- anole_adm_cv
  anole_adm_no_group$group <- NULL
  expect_error(gt_admix_reorder_q(anole_adm_no_group),
               "You must provide a group variable if there is no grouping information in the gt_admix object")
  # use gt_adm without id and group
  anole_adm_no_group$id <- NULL
  anole_adm_no_group$group <- NULL
  adm_reord_no_meta <- gt_admix_reorder_q(anole_adm_no_group, group=anole_gt$population)
  # now check that we have an id and group elements in the list
  expect_true(length(adm_reord_no_meta$id)==nrow(anole_gt))
  expect_true(length(adm_reord_no_meta$group)==nrow(anole_gt))
})
