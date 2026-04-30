bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
missing_gt <- gen_tibble(
  bed_file,
  backingfile = tempfile("missing_"),
  quiet = TRUE
)
missing_gt <- gt_impute_simple(missing_gt, method = "mode")
missing_pca <- missing_gt %>% gt_pca_partialSVD()

test_that("gt_pcadapt works on gt_pca object", {
  missing_pcadapt <- gt_pcadapt(missing_gt, missing_pca, k = 3)
  expect_true((inherits(missing_pcadapt, "gt_pcadapt")))
})

test_that("n_cores can be set", {
  expect_true(getOption("bigstatsr.check.parallel.blas"))
  one_core <- gt_pcadapt(missing_gt, missing_pca, k = 3, n_cores = 1)
  expect_true(getOption("bigstatsr.check.parallel.blas"))
  two_core <- gt_pcadapt(missing_gt, missing_pca, k = 3, n_cores = 2)
  expect_true(getOption("bigstatsr.check.parallel.blas"))
})

test_that("gt_pcadapt on gen_tbl vcf",{
  vcf_path <-
    system.file("/extdata/anolis/punctatus_t70_s10_n46_filtered.recode.vcf.gz",
                package = "tidypopgen"
    )
  anolis <- gen_tibble(
    vcf_path,
    quiet = TRUE,
    backingfile = tempfile(),
    parser = "vcfR"
  )
  anolis <- gt_impute_simple(anolis, method = "mode")
  anolis <- anolis %>% select_loci_if(loci_maf(genotypes) > 0)
  anolis <- gt_update_backingfile(anolis)

  dim(.gt_get_fbm(anolis))
  gt_uses_imputed(anolis)
  gt_has_imputed(anolis)
  gt_set_imputed(anolis, TRUE)

  pop_b_vcf_pca <- gt_pca_randomSVD(anolis)
  pop_b_vcf_pcadapt <- gt_pcadapt(anolis, pop_b_vcf_pca, k = 2)
})
