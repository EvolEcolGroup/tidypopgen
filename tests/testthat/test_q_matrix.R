test_that("q_matrix reads q_files", {
  # read a single .Q
  q_path <- system.file(
    "/extdata/anolis/anolis_ld_run1.3.Q",
    package = "tidypopgen"
  )

  # read from df
  q_mat <- read.table(q_path)
  anolis_q_k3_mat <- q_matrix(q_mat)
  expect_true(inherits(anolis_q_k3_mat, "q_matrix"))

  # read from matrix
  q_mat <- as.matrix(q_mat)
  anolis_q_k3_mat <- q_matrix(q_mat)
  expect_true(inherits(anolis_q_k3_mat, "q_matrix"))

  # read multiple .Q
  q_folder <- system.file("/extdata/anolis", package = "tidypopgen")
  anolis_q <- read_q_files(q_folder)
  expect_true(inherits(anolis_q, "gt_admix"))
  expect_true(inherits(anolis_q$Q[[1]], "q_matrix"))

  # check errors
  non_path <- "an/invalid/path"
  expect_error(read_q_files(non_path), "Input is not a valid directory")
  path_no_q <- system.file("/extdata/", package = "tidypopgen")
  expect_error(read_q_files(path_no_q), "No .Q files found in the directory")
})

test_that("get_q_matrix returns correct q-matrix", {
  # read multiple q into a list
  q_folder <- system.file("/extdata/anolis", package = "tidypopgen")
  q_list <- read_q_files(q_folder)

  # check returns a single q-matrix object
  anolis_q_k3_r1 <- get_q_matrix(q_list, k = 3, run = 1)
  expect_true(inherits(anolis_q_k3_r1, "q_matrix"))

  # check number of cols of q-matrix is equal to k
  expect_equal(ncol(anolis_q_k3_r1), 3)
  anolis_q_k4_r1 <- get_q_matrix(q_list, k = 4, run = 1)
  expect_equal(ncol(anolis_q_k4_r1), 4)

  # check errors if outside of k or run range
  expect_error(
    get_q_matrix(q_list, k = 5, run = 1),
    "Specified k value not found in gt_admix"
  )

  # check errors if give non 'gt_admix' object
  q_list_no_structure <- list(
    as.matrix(read.table(
      system.file("/extdata/anolis/anolis_ld_run1.3.Q", package = "tidypopgen"),
      header = TRUE
    )),
    as.matrix(read.table(
      system.file("/extdata/anolis/anolis_ld_run1.4.Q", package = "tidypopgen"),
      header = TRUE
    ))
  )
  expect_error(
    get_q_matrix(q_list_no_structure, k = 3, run = 1),
    "Input must be a 'gt_admix' object"
  )
  expect_error(
    get_q_matrix(q_list$Q[[1]], k = 3, run = 1),
    "Input must be a 'gt_admix' object"
  )
  expect_error(
    get_q_matrix(q_list, k = 3, run = 3),
    "Specified run is out of range"
  )
})

test_that("get_p_matrix returns correct p-matrix", {
  # read multiple q into a list
  q_folder <- system.file("/extdata/anolis", package = "tidypopgen")
  q_list <- read_q_files(q_folder)

  # check error if 'q_matrix_list' object doesn't contain p-files
  expect_error(
    get_p_matrix(q_list, k = 5, run = 1),
    "Input object does not contain any P matrices"
  )

  # check errors if give non 'q_matrix_list' object
  q_list_no_structure <- list(
    as.matrix(read.table(
      system.file("/extdata/anolis/anolis_ld_run1.3.Q", package = "tidypopgen"),
      header = TRUE
    )),
    as.matrix(read.table(
      system.file("/extdata/anolis/anolis_ld_run1.4.Q", package = "tidypopgen"),
      header = TRUE
    ))
  )
  expect_error(
    get_p_matrix(q_list_no_structure, k = 3, run = 1),
    "Input must be a 'gt_admix' object"
  )
  expect_error(
    get_p_matrix(q_list$Q[[1]], k = 3, run = 1),
    "Input must be a 'gt_admix' object"
  )

  # add P files to q_list
  q_list$P <- list(
    as.matrix(read.table(
      system.file("/extdata/anolis/anolis_ld_run1.3.P", package = "tidypopgen"),
      header = TRUE
    )),
    as.matrix(read.table(
      system.file("/extdata/anolis/anolis_ld_run1.4.P", package = "tidypopgen"),
      header = TRUE
    ))
  )

  # check returns a single q-matrix object
  anolis_p_k3_r1 <- get_p_matrix(q_list, k = 3, run = 1)
  expect_true(inherits(anolis_p_k3_r1, "matrix"))

  # check number of cols of p-matrix is equal to k
  anolis_p_k3_r1 <- get_p_matrix(q_list, k = 3, run = 1)
  anolis_p_k4_r1 <- get_p_matrix(q_list, k = 4, run = 1)
  expect_equal(ncol(anolis_p_k3_r1), 3)
  expect_equal(ncol(anolis_p_k4_r1), 4)

  # check errors if outside of k or run range
  expect_error(
    get_p_matrix(q_list, k = 5, run = 1),
    "Specified k value not found in gt_admix"
  )
  expect_error(
    get_p_matrix(q_list, k = 3, run = 3),
    "Specified run is out of range"
  )
})

if (rlang::is_installed("readr")) {
  test_that("tidying and augmenting a q_matrix", {
    # read a single .Q
    q_path <- system.file(
      "/extdata/anolis/anolis_ld_run1.3.Q",
      package = "tidypopgen"
    )
    q_mat <- read.table(q_path)
    anolis_q_k3_mat <- q_matrix(q_mat)

    # create gt and metadata
    vcf_path <- system.file(
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
    anole_gt <- anole_gt %>%
      mutate(population = pops$pop[match(pops$ID, .data$id)])

    # tidy without group info
    tidy_q <- tidy(anolis_q_k3_mat, anole_gt)
    expect_false("group" %in% colnames(tidy_q))
    expect_equal(tidy_q$id, rep(anole_gt$id, each = 3))

    # augment without group info
    augment_q <- augment(anolis_q_k3_mat, anole_gt)
    expect_false("group" %in% colnames(anolis_q_k3_mat))
    expect_equal(augment_q$id, anole_gt$id)
    expect_true(inherits(augment_q, "gen_tbl"))

    anole_gt <- anole_gt %>% group_by(population)

    # tidy using group from gen_tibble
    tidy_q <- tidy(anolis_q_k3_mat, anole_gt)
    expect_true("group" %in% colnames(tidy_q))
    expect_equal(tidy_q$group, rep(anole_gt$population, each = 3))
    expect_equal(tidy_q$id, rep(anole_gt$id, each = 3))

    # augment using group from gen_tibble
    augment_q <- augment(anolis_q_k3_mat, anole_gt)
    expect_true(inherits(augment_q, "gen_tbl"))
    expect_true("group" %in% colnames(augment_q))
    expect_equal(augment_q$group, anole_gt$population)
    expect_equal(augment_q$id, anole_gt$id)
    # the grouped class should be preserved
    expect_equal(class(augment_q), class(anole_gt))
    expect_equal(group_vars(augment_q), group_vars(anole_gt))
    # and we should be able to filter the augmented gen_tibble
    pop_to_keep <- unique(anole_gt$population)[1]
    augment_sub <- augment_q %>% filter(population == pop_to_keep)
    expect_equal(class(augment_sub), class(anole_gt))
    expect_true(all(augment_sub$population == pop_to_keep))
  })
}

test_that("tidy and augment use the id attribute of a q_matrix", {
  example_gt <- load_example_gt("gen_tbl")
  n_ind <- nrow(example_gt)
  q <- seq(0.1, 0.9, length.out = n_ind)
  admix_obj <- structure(
    list(
      k = 2L,
      Q = list(q_matrix(matrix(c(q, 1 - q), ncol = 2))),
      id = example_gt$id,
      algorithm = "ADMIXTURE"
    ),
    class = c("gt_admix", "list")
  )
  # a q_matrix with an id attribute, but no group attribute
  q_mat <- get_q_matrix(admix_obj, k = 2, run = 1)
  expect_null(attr(q_mat, "group"))
  # augment with data in a different order
  example_gt_rev <- example_gt %>% dplyr::arrange(dplyr::desc(id))
  augment_q <- augment(q_mat, data = example_gt_rev)
  expect_equal(augment_q$id, example_gt_rev$id)
  expect_identical(
    augment_q$.Q1,
    q[match(example_gt_rev$id, example_gt$id)]
  )
  # tidy with data in a different order
  tidy_q <- tidy(q_mat, data = example_gt_rev)
  expect_equal(tidy_q$id, rep(example_gt$id, each = 2))
  expect_equal(tidy_q$percentage[tidy_q$q == ".Q1"], q)
})
