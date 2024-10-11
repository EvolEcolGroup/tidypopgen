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
  package_path <- find.package("tidypopgen")
  folder <- "/extdata/anolis"
  q_folder <- file.path(package_path, folder)
  anolis_q <- q_matrix(q_folder)
  expect_true(inherits(anolis_q,"q_matrix_list"))
  expect_true(inherits(anolis_q[[1]][[1]],"q_matrix"))

  #check errors
  P_path <- system.file("/extdata/anolis/anolis_ld_run1.3.P", package = "tidypopgen")
  expect_error(q_matrix(P_path),"Input file does not end in '.Q'")
  non_path <- "an/invalid/path"
  expect_error(q_matrix(non_path),"Input is not a valid dataframe, file or directory")
  path_no_q <- system.file("/extdata/", package = "tidypopgen")
  expect_error(q_matrix(path_no_q),"No .Q files found in the directory")
})

