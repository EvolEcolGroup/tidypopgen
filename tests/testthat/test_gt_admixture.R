test_that("run admixture correctly", {
  # skip if admixture is not installed
  skip_if((system2("which", args = "admixture") != 0)||!requireNamespace("fastmixturer", quietly = TRUE))
  # set the input file
  vcf_path <- system.file("/extdata/anolis/punctatus_t70_s10_n46_filtered.recode.vcf.gz",
                          package = "tidypopgen")
  anole_gt <- gen_tibble(vcf_path, quiet = TRUE, backingfile = tempfile("anolis_"))
  pops_path <- system.file("/extdata/anolis/plot_order_punctatus_n46.csv",
                           package = "tidypopgen")
  pops <- readr::read_csv(pops_path, show_col_types = FALSE)
  anole_gt <- anole_gt %>% mutate(id = gsub('punc_',"",.data$id,))
  anole_gt <- anole_gt %>% mutate(population = pops$pop[match(pops$ID,.data$id)])
 # anole_gt <- anole_gt %>% group_by(population)
  # we create a plink file to test the function
  anole_plink <- gt_as_plink(anole_gt, file = tempfile(), chromosomes_as_int=TRUE)
  # run admixture
  gt_admixture(anole_plink, k = 3, crossval = FALSE, n_cores = 1, conda_env = "none")

})
