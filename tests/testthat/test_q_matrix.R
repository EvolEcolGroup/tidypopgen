test_that("q_matrix reads q_files",{
  #read a single .Q
  Q_path <- system.file("/extdata/anolis/anolis_ld_run1.3.Q", package = "tidypopgen")
  anolis_q_k3 <- q_matrix(Q_path)
  expect_true(inherits(anolis_q_k3,"q_matrix"))

  #read from df
  q_mat <- read.table(Q_path)
  anolis_q_k3_mat <- q_matrix(q_mat)
  expect_true(inherits(anolis_q_k3_mat,"q_matrix"))

  #read from matrix
  q_mat <- as.matrix(q_mat)
  anolis_q_k3_mat <- q_matrix(q_mat)
  expect_true(inherits(anolis_q_k3_mat,"q_matrix"))

  #read multiple .Q
  q_folder <- system.file("/extdata/anolis", package = "tidypopgen")
  anolis_q <- q_matrix(q_folder)
  expect_true(inherits(anolis_q,"q_matrix_list"))
  expect_true(inherits(anolis_q$q_matrices[[1]],"q_matrix"))

  #check errors
  P_path <- system.file("/extdata/anolis/anolis_ld_run1.3.P", package = "tidypopgen")
  expect_error(q_matrix(P_path),"Input file does not end in '.Q'")
  non_path <- "an/invalid/path"
  expect_error(q_matrix(non_path),"Input is not a valid dataframe, file, directory, or list of q-matrices")
  path_no_q <- system.file("/extdata/", package = "tidypopgen")
  expect_error(q_matrix(path_no_q),"No .Q files found in the directory")
})

test_that("get_q returns correct matrix",{
  #read multiple q into a list
  q_folder <- system.file("/extdata/anolis", package = "tidypopgen")
  q_list <- q_matrix(q_folder)

  #check returns a single q-matrix object
  anolis_q_k3_r1 <- get_q(q_list, k=3, run=1)
  expect_true(inherits(anolis_q_k3_r1,"q_matrix"))

  #check number of cols of q-matrix is equal to k
  expect_equal(ncol(anolis_q_k3_r1),3)
  anolis_q_k4_r1 <- get_q(q_list, k=4, run=1)
  expect_equal(ncol(anolis_q_k4_r1),4)

  #check errors if outside of k or run range
  expect_error(get_q(q_list, k=5, run=1),"Specified k value not found in q_matrix_list")
  expect_error(get_q(q_list, k=3, run=2),"Specified run is out of range for the given k value")

  #check errors if give non 'q_matrix_list' object
  q_list_no_structure <- list(as.matrix(read.table(system.file("/extdata/anolis/anolis_ld_run1.3.Q", package = "tidypopgen"), header = TRUE)),
                              as.matrix(read.table(system.file("/extdata/anolis/anolis_ld_run1.4.Q", package = "tidypopgen"), header = TRUE)))
  expect_error(get_q(q_list_no_structure, k=3, run=1),"Input must be a 'q_matrix_list' object")
  expect_error(get_q(q_list$q_matrices[[1]], k=3, run=1),"Input must be a 'q_matrix_list' object")
})

test_that("get_p only returns p",{
  #read multiple q into a list
  q_folder <- system.file("/extdata/anolis", package = "tidypopgen")
  q_list <- q_matrix(q_folder)

  #check error if 'q_matrix_list' object doesn't contain p-files
  expect_error(get_p(q_list, k=5, run=1),"Input object does not contain any P matrices")

  #check errors if give non 'q_matrix_list' object
  q_list_no_structure <- list(as.matrix(read.table(system.file("/extdata/anolis/anolis_ld_run1.3.Q", package = "tidypopgen"), header = TRUE)),
                              as.matrix(read.table(system.file("/extdata/anolis/anolis_ld_run1.4.Q", package = "tidypopgen"), header = TRUE)))
  expect_error(get_p(q_list_no_structure, k=3, run=1),"Input must be a 'q_matrix_list' object")
  expect_error(get_p(q_list$q_matrices[[1]], k=3, run=1),"Input must be a 'q_matrix_list' object")
})

test_that("summary of q_matrix",{
  #read multiple q into a list
  q_folder <- system.file("/extdata/anolis", package = "tidypopgen")
  q_list <- q_matrix(q_folder)
  sum <- summary.q_matrix_list(q_list)

  #make expected
  expected_summary <- data.frame(
    K = c(3, 4),
    Repeats = c(1, 1))

  #check summary of q_matrix_list
  expect_equal(sum, expected_summary)

  #check errors if give non 'q_matrix_list' object
  q_list_no_structure <- list(as.matrix(read.table(system.file("/extdata/anolis/anolis_ld_run1.3.Q", package = "tidypopgen"), header = TRUE)),
                              as.matrix(read.table(system.file("/extdata/anolis/anolis_ld_run1.4.Q", package = "tidypopgen"), header = TRUE)))
  expect_error(summary.q_matrix_list(q_list_no_structure),"Input must be a 'q_matrix_list' object")

})

