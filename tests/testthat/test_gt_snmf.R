# limit number of threads for tests
data.table::setDTthreads(2)
if (rlang::is_installed("RhpcBLASctl")) {
  RhpcBLASctl::blas_set_num_threads(2)
  RhpcBLASctl::omp_set_num_threads(2)
}

skip_if_not_installed("LEA")

# set the input file
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
anole_gt <- anole_gt %>% mutate(population = pops$pop[match(pops$ID, .data$id)])

library(LEA)

test_that("gt_snmf error messages", {
  # wrong file input
  expect_error(
    gt_snmf(
      x = pops_path,
      project = "force",
      k = 1:10,
      entropy = TRUE,
      percentage = 0.5,
      n_runs = 1,
      seed = 1,
      alpha = 100
    ),
    "The input file must be a .geno"
  )
  # file that doesn't exist
  invalid_geno <- paste0(tempfile(), ".geno")
  expect_error(
    gt_snmf(
      x = invalid_geno,
      project = "force",
      k = 1:10,
      entropy = TRUE,
      percentage = 0.5,
      n_runs = 1,
      seed = 1,
      alpha = 100
    ),
    "does not exist"
  )
  # random object
  expect_error(
    gt_snmf(
      x = 1,
      project = "force",
      k = 1:10,
      entropy = TRUE,
      percentage = 0.5,
      n_runs = 1,
      seed = 1,
      alpha = 100
    ),
    "x must be a gen_tibble or a character giving the path"
  )
  expect_error(
    gt_snmf(
      x = anole_gt,
      project = "force",
      k = 1:10,
      entropy = TRUE,
      percentage = 0.5,
      n_runs = 1,
      seed = c(1, 2),
      alpha = 100
    ),
    "'seed' should be a vector of length 'n_runs'"
  )
})

test_that("gt_snmf from file and from gen_tibble are the same", {
  # using.geno file
  input_file <- gt_as_geno_lea(anole_gt)
  capture.output(anole_snmf_file <- gt_snmf(
    x = input_file,
    project = "force",
    k = 1:10,
    entropy = TRUE,
    percentage = 0.5,
    n_runs = 1,
    seed = 1,
    alpha = 100
  ))
  # using gen_tibble
  capture.output(anole_snmf_gt <- gt_snmf(
    x = anole_gt,
    project = "force",
    k = 1:10,
    entropy = TRUE,
    percentage = 0.5,
    n_runs = 1,
    seed = 1,
    alpha = 100
  ))
  # check that the results are the same
  expect_equal(anole_snmf_file$Q, anole_snmf_gt$Q)
  expect_equal(anole_snmf_file$G, anole_snmf_gt$G)
  expect_equal(anole_snmf_file$cv, anole_snmf_gt$cv)

  # check against snmf
  capture.output(anole_snmf_lea <- LEA::snmf(
    input.file = input_file,
    project = "force",
    K = 1:10,
    entropy = TRUE,
    percentage = 0.5,
    repetitions = 1,
    seed = 1,
    alpha = 100
  ))
  lea_q <- LEA::Q(anole_snmf_lea, K = 1, run = 1)
  gt_q <- get_q_matrix(anole_snmf_file, k = 1, run = 1)
  expect_true(all(lea_q) == all(gt_q))

  # after removing entropy arguments
  anole_snmf_gt_ne <- gt_snmf(
    x = anole_gt,
    project = "force",
    k = 1:10,
    n_runs = 1,
    seed = 1,
    alpha = 100
  )
  # the results differ from the runs with entropy (which are estimated on
  # the masked genotypes), so we check that we get the latest run (run 4,
  # after the 3 runs above) from the project
  expect_equal(anole_snmf_gt_ne$k, anole_snmf_gt$k)
  expect_null(anole_snmf_gt_ne$cv)
  anole_project <- LEA::load.snmfProject(
    sub("\\.geno$", ".snmfProject", input_file)
  )
  expect_equal(
    unclass(anole_snmf_gt_ne$Q[[3]]),
    LEA::Q(anole_project, K = 3, run = 4),
    ignore_attr = TRUE
  )
  expect_equal(
    unclass(anole_snmf_gt_ne$G[[3]]),
    LEA::G(anole_project, K = 3, run = 4),
    ignore_attr = TRUE
  )

  # test plot
  plt <- autoplot(anole_snmf_gt)
  expect_true(plt$labels$y == "Cross-Entropy")
})

test_that("gt_snmf cross-entropy matches k and run", {
  example_gt <- load_example_gt("gen_tbl")
  input_file <- gt_as_geno_lea(example_gt)
  capture.output(example_snmf <- gt_snmf(
    x = input_file,
    k = 2:3,
    project = "new",
    n_runs = 2,
    seed = c(1, 2),
    entropy = TRUE,
    percentage = 0.5,
    alpha = 100
  ))
  snmf_project <- LEA::load.snmfProject(
    sub("\\.geno$", ".snmfProject", input_file)
  )
  expect_equal(example_snmf$k, c(2, 2, 3, 3))
  lea_cv <- c(
    LEA::cross.entropy(snmf_project, K = 2),
    LEA::cross.entropy(snmf_project, K = 3)
  )
  expect_equal(example_snmf$cv, lea_cv, tolerance = 1e-5)
})

test_that("gt_snmf default seed works with several runs", {
  example_gt <- load_example_gt("gen_tbl")
  capture.output(example_snmf <- gt_snmf(
    x = example_gt,
    k = 2,
    project = "new",
    n_runs = 2,
    alpha = 100
  ))
  expect_equal(example_snmf$k, c(2, 2))
  expect_length(example_snmf$Q, 2)
  # a single seed is still not accepted for several runs
  expect_error(
    gt_snmf(
      x = example_gt,
      k = 2,
      project = "new",
      n_runs = 2,
      seed = 1,
      alpha = 100
    ),
    "'seed' should be a vector of length 'n_runs'"
  )
})

test_that("gt_snmf returns the runs of the current call in a project", {
  example_gt <- load_example_gt("gen_tbl")
  input_file <- gt_as_geno_lea(example_gt)
  project_file <- sub("\\.geno$", ".snmfProject", input_file)
  capture.output(snmf_1 <- gt_snmf(
    x = input_file,
    k = 2,
    project = "new",
    n_runs = 1,
    seed = 11,
    alpha = 100
  ))
  # add runs to the existing project
  capture.output(snmf_2 <- gt_snmf(
    x = input_file,
    k = 2:3,
    project = "continue",
    n_runs = 2,
    seed = c(5, 6),
    entropy = TRUE,
    percentage = 0.5,
    alpha = 1
  ))
  snmf_project <- LEA::load.snmfProject(project_file)
  expect_equal(
    unclass(snmf_1$Q[[1]]),
    LEA::Q(snmf_project, K = 2, run = 1),
    ignore_attr = TRUE
  )
  # the first call created run 1 of K = 2, so the second call created runs 2
  # and 3 for K = 2, and runs 1 and 2 for K = 3
  lea_runs <- list(c(2, 2), c(2, 3), c(3, 1), c(3, 2))
  for (i in seq_along(lea_runs)) {
    expect_equal(
      unclass(snmf_2$Q[[i]]),
      LEA::Q(snmf_project, K = lea_runs[[i]][1], run = lea_runs[[i]][2]),
      ignore_attr = TRUE
    )
    expect_equal(
      unclass(snmf_2$G[[i]]),
      LEA::G(snmf_project, K = lea_runs[[i]][1], run = lea_runs[[i]][2]),
      ignore_attr = TRUE
    )
    expect_equal(
      snmf_2$cv[i],
      LEA::cross.entropy(snmf_project,
        K = lea_runs[[i]][1],
        run = lea_runs[[i]][2]
      )[1, 1],
      tolerance = 1e-5
    )
  }
})
