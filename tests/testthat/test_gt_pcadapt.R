# limit number of threads for tests
data.table::setDTthreads(2)
if (rlang::is_installed("RhpcBLASctl")) {
  RhpcBLASctl::blas_set_num_threads(2)
  RhpcBLASctl::omp_set_num_threads(2)
}

bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
missing_gt <- gen_tibble(
  bed_file,
  backingfile = tempfile("missing_"),
  quiet = TRUE
)
missing_gt <- gt_impute_simple(missing_gt, method = "mode")
# the PCA of the imputed genotypes
gt_set_imputed(missing_gt, set = TRUE)
missing_pca <- missing_gt %>% gt_pca(n_cores = 2)
gt_set_imputed(missing_gt, set = FALSE)

test_that("gt_pcadapt works on gt_pca object", {
  missing_pcadapt <- gt_pcadapt(missing_gt, missing_pca, k = 3)
  expect_true((inherits(missing_pcadapt, "gt_pcadapt")))
})

test_that("gt_pcadapt with a PCA of missing genotypes", {
  # a PCA that compared individuals over the loci genotyped in both
  pairwise_pca <- missing_gt %>% gt_pca(k = 3, n_cores = 2)
  expect_true(attr(pairwise_pca, "missing"))
  pairwise_pcadapt <- gt_pcadapt(missing_gt, pairwise_pca, k = 1)
  # each locus is tested on the individuals genotyped at it (not on the
  # imputed genotypes of missing_gt): with k = 1, the statistic is the
  # squared deviation from the median of the t-score of the regression of the
  # genotypes on the scores
  expect_false(gt_uses_imputed(missing_gt))
  geno <- show_genotypes(missing_gt)
  u <- pairwise_pca$u[, 1]
  t_scores <- apply(geno, 2, function(g) {
    obs <- !is.na(g)
    r <- stats::cor(g[obs], u[obs])
    r * sqrt((sum(obs) - 2) / (1 - r^2))
  })
  expect_equal(
    pairwise_pcadapt$score,
    (t_scores - stats::median(t_scores))^2
  )
  # loci with the same genotype in every individual genotyped are named
  het_gt <- gen_tibble(
    x = rbind(
      c(2, 1, 0, 1, 1, 1), c(0, 1, 2, 0, 1, NA), c(1, 2, 0, 1, 0, 1),
      c(2, 0, 1, 2, 1, 1), c(0, 2, 1, 0, 2, NA), c(1, 1, 2, 2, 0, 1),
      c(2, 0, 0, 1, 2, 1), c(0, 1, 1, 0, 1, 1)
    ),
    indiv_meta = data.frame(
      id = paste0("ind", 1:8), population = rep(c("a", "b"), each = 4)
    ),
    loci = data.frame(
      name = paste0("rs", 1:6), chromosome = 1, position = 1:6 * 100,
      genetic_dist = 0, allele_ref = "A", allele_alt = "T"
    ),
    backingfile = tempfile(), quiet = TRUE
  )
  het_pca <- gt_pca(het_gt, k = 2, n_cores = 1)
  expect_error(
    gt_pcadapt(het_gt, het_pca, k = 2),
    "1 loci are monomorphic or heterozygous in every individual \\(genotyped"
  )
})

test_that("gt_pcadapt checks that individuals match the pca", {
  # same individuals, in a different order
  rev_gt <- missing_gt[rev(seq_len(nrow(missing_gt))), ]
  expect_error(
    gt_pcadapt(rev_gt, missing_pca, k = 3),
    "the individuals in 'x' do not correspond to those in 'pca'"
  )
})

test_that("gt_pcadapt looks for missing genotypes only where they are used", {
  # the 133 complete loci of example-missing, without resaving the backing
  # matrix (where the other loci keep their missing genotypes)
  raw_gt <- gen_tibble(
    bed_file,
    backingfile = tempfile("missing_"),
    quiet = TRUE
  )
  complete_gt <- raw_gt %>% select_loci_if(loci_missingness(genotypes) == 0)
  complete_pca <- gt_pca(complete_gt, k = 3, n_cores = 1)
  # bigsnpr's check sampled the whole backing matrix, and stopped with
  # "You can't have missing values"
  set.seed(1)
  expect_s3_class(gt_pcadapt(complete_gt, complete_pca, k = 3), "gt_pcadapt")
  # missing genotypes in the loci used are still found
  expect_error(
    gt_pcadapt(raw_gt, complete_pca, k = 3),
    "gt_pcadapt\\(\\) needs complete or imputed genotypes"
  )
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
  anolis_pca <- gt_pca(anolis, method = "random", n_cores = 1)
  # fails in gt_pcadapt (snp_pcadapt() reported "You can't have missing
  # values in 'U'"): the loci heterozygous in every individual are found
  expect_error(
    gt_pcadapt(anolis, anolis_pca, k = 2),
    "27 loci are monomorphic or heterozygous in every individual"
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
  x <- pcadapt::pcadapt(input = filename, K = 2)
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
  anolis_pca <- gt_pca(anolis, method = "random", n_cores = 1)
  # and now it works
  anolis_pcadapt <- gt_pcadapt(anolis, anolis_pca, k = 2)
  expect_true(inherits(anolis_pcadapt, "gt_pcadapt"))
})
