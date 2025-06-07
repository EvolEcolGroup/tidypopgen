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
  expect_equal(anole_snmf_gt$Q, anole_snmf_gt_ne$Q)
  expect_equal(anole_snmf_gt$G, anole_snmf_gt_ne$G)

  # test plot
  plt <- autoplot(anole_snmf_gt)
  expect_true(plt$labels$y == "Cross-Entropy")
})
