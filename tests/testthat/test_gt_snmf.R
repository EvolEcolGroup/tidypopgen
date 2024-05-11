vcf_path <- system.file("/extdata/anolis/punctatus_t70_s10_n46_filtered.recode.vcf.gz",
                        package = "tidypopgen")
anole_gt <- gen_tibble(vcf_path, quiet = TRUE, backingfile = tempfile())
test_that("snmf clusters correctly",{
  geno_file <- gt_as_geno_lea(anole_gt)
  # silence verbose snmf function
  sink(file=tempfile())
  snmf_project <- LEA::snmf(geno_file, K = 3)
  # end of sinking
  sink()
  expect_true(inherits(snmf_project,"snmfProject"))
})
