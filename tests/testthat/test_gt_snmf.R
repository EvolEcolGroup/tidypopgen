vcf_path <- system.file("/extdata/anolis/punctatus_t70_s10_n46_filtered.recode.vcf.gz",
                        package = "tidypopgen")
bed_path <- gt_vcf_to_bed(vcf_path, bed_path = tempfile("anolis_"))
anole_gt <- gen_tibble(bed_path, quiet = TRUE)
test_that("snmf clusters correctly",{
  geno_file <- gt_write_lea_geno(anole_gt)
  # silence verbose snmf function
  sink(file=tempfile())
  snmf_project <- LEA::snmf(geno_file, K = 3)
  # end of sinking
  sink()
  expect_true(inherits(snmf_project,"snmfProject"))
})
