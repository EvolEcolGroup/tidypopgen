test_that("run admixture correctly", {
  # skip if admixture is not installed
  skip_if((system2("which", args = "admixture", stdout = NULL) != 0)||!requireNamespace("fastmixturer", quietly = TRUE))
  # set the input file
  vcf_path <- system.file("/extdata/anolis/punctatus_t70_s10_n46_filtered.recode.vcf.gz",
                          package = "tidypopgen")
  anole_gt <- gen_tibble(vcf_path, quiet = TRUE, backingfile = tempfile("anolis_"))
  pops_path <- system.file("/extdata/anolis/plot_order_punctatus_n46.csv",
                           package = "tidypopgen")
  pops <- readr::read_csv(pops_path, show_col_types = FALSE)
  anole_gt <- anole_gt %>% mutate(id = gsub('punc_',"",.data$id,))
  anole_gt <- anole_gt %>% mutate(population = pops$pop[match(pops$ID,.data$id)])
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
  anole_adm3 <- gt_admixture(anole_plink, k = 3, crossval = TRUE, n_cores = 2, seed = 345, conda_env = "none")
  expect_false(is.null(anole_adm3$cv))
  anole_adm4 <- gt_admixture(anole_plink, k = 3, crossval = TRUE, n_cores = 2, seed = 123, conda_env = "none")
  anole_adm_comb2 <- c(anole_adm3, anole_adm4)


})
