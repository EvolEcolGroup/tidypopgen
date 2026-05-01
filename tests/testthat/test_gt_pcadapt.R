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

skip_if_not_installed("pcadapt")

test_that("gt_pcadapt 'U' error with heavily imputed SNPs", {
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
  # remove monomorphic SNPs
  anolis <- anolis %>% select_loci_if(loci_maf(genotypes) > 0.3)
  # update backingfile
  anolis <- gt_update_backingfile(anolis, quiet = TRUE)
  # impute
  anolis <- gt_impute_simple(anolis, method = "mode")
  # set imputed to TRUE (just to be sure)
  gt_set_imputed(anolis, TRUE)
  # test pca
  anolis_pca <- gt_pca_randomSVD(anolis)
  # fails in gt_pcadapt
  expect_error(
    gt_pcadapt(anolis, anolis_pca, k = 2),
    "You can't have missing values in 'U'"
  )

  # double check our data is imputed (no missingness)
  expect_true(gt_has_imputed(anolis))
  expect_true(inherits(attributes(anolis$genotypes)$fbm, "FBM.code256"))
  # double check the pca scores do not contain NA's
  expect_false(any(is.na(anolis_pca$u)))

  # Now we write the same data out as a bed file
  path <- file.path(tempdir(), "anolis_plink")
  gt_as_plink(anolis, file = path)
  library(pcadapt)
  filename <- read.pcadapt(paste0(path, ".bed"), type = "bed")
  # and pcadapt works on this bed file
  # (although, this uses pcadapt, while our function uses snp_pcadapt)
  x <- pcadapt(input = filename, K = 2)
  expect_true(inherits(x, "pcadapt"))

  # next, we used bigsnpr:::multiLinReg directly to generate the tscores object
  # this shows that NA's are introduced
  # we find NA-causing SNPs from tscores
  # [1]   36   85  182  184  197  265  381  411  471  472  473  530  576  706
  # 713  762  823  846  854  855  856
  # [22]  873  878  886 1002 1073 1081
  # and write these as a vector to remove
  snps_to_remove <- c(
    36, 85, 182, 184, 197, 265, 381, 411, 471, 472, 473, 530, 576, 706,
    713, 762, 823, 846, 854, 855, 856,
    873, 878, 886, 1002, 1073, 1081
  )
  anolis_bad_snps <- anolis %>% select_loci(all_of(snps_to_remove))
  # all of these SNPs have a MAF of 0.5 (imputed)
  loci_maf(anolis_bad_snps)
  gt_set_imputed(anolis_bad_snps, TRUE)
  show_genotypes(anolis_bad_snps)

  # remove these SNPs
  anolis <- anolis %>% select_loci(-all_of(snps_to_remove))
  # check MAF is sensible
  loci_maf(anolis)
  # update backingfile
  anolis <- gt_update_backingfile(anolis, quiet = TRUE)
  # test pca
  anolis_pca <- gt_pca_randomSVD(anolis)
  # and now it works
  anolis_pcadapt <- gt_pcadapt(anolis, anolis_pca, k = 2)
  expect_true(inherits(anolis_pcadapt, "gt_pcadapt"))
})
