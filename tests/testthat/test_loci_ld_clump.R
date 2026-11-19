test_indiv_meta <- data.frame(
  id = c("a", "b", "c"),
  population = c("pop1", "pop1", "pop2")
)
test_genotypes <- rbind(
  c(2, 2, 0, 1, 1, 2),
  c(2, 2, 0, 0, 0, 1),
  c(2, 2, 0, 0, 1, 1)
)
test_loci <- data.frame(
  name = paste0("rs", 1:6),
  chromosome = c(1, 1, 1, 1, 2, 2),
  position = c(3, 5, 65, 343, 23, 456),
  genetic_dist = as.double(rep(0, 6)),
  allele_ref = c("A", "T", "C", "G", "C", "T"),
  allele_alt = c("T", "C", NA, "C", "G", "A")
)


test_gt <- gen_tibble(
  x = test_genotypes,
  indiv_meta = test_indiv_meta,
  loci = test_loci,
  backingfile = tempfile(),
  quiet = TRUE
)

test_that("ld clumping runs", {
  keep <- loci_ld_clump(test_gt, thr_r2 = 0.2, return_id = TRUE)
  expect_true(all.equal(keep, c(1, 2, 3, 4, 6)) == TRUE)
})

test_that("loci_ld_clump with missing genotypes", {
  # the PLINK fixture: 60 individuals, two chromosomes of 100 loci in LD, with
  # missing genotypes (see data-raw/validate/create_plink_ld_test_data.R)
  ld_gt <- gen_tibble(test_path("testdata", "plink_ld", "ld.bed"),
    backingfile = tempfile("ld_"), quiet = TRUE
  )
  geno <- show_genotypes(ld_gt)
  loci <- show_loci(ld_gt)
  # the reference: bigsnpr's algorithm, sequentially, with r2 over the
  # individuals genotyped at both loci (loci in order of decreasing MAF; a
  # locus is kept unless linked to a kept locus within 50 kb)
  maf <- loci_maf(ld_gt)
  ref <- rep(FALSE, ncol(geno))
  for (chrom in unique(loci$chromosome)) {
    idx <- which(loci$chromosome == chrom)
    kept <- integer(0)
    for (j in idx[order(maf[idx], decreasing = TRUE)]) {
      linked <- FALSE
      for (k in kept[abs(loci$position[kept] - loci$position[j]) <= 50000]) {
        shared <- !is.na(geno[, j]) & !is.na(geno[, k])
        r2 <- suppressWarnings(stats::cor(geno[shared, j], geno[shared, k])^2)
        if (sum(shared) >= 2 && !is.na(r2) && r2 > 0.2) {
          linked <- TRUE
          break
        }
      }
      if (!linked) kept <- c(kept, j)
    }
    ref[kept] <- TRUE
  }
  expect_identical(loci_ld_clump(ld_gt, thr_r2 = 0.2, size = 50), ref)
  expect_identical(
    loci_ld_clump(ld_gt, thr_r2 = 0.2, size = 50, n_cores = 2),
    ref
  )
  # with a minimum above the number of individuals, no pair has an r2: with
  # "keep", every locus is kept; with "remove", every pair is linked, so that
  # no two loci kept on a chromosome are within 50 kb
  keep_all <- loci_ld_clump(ld_gt, thr_r2 = 0.2, size = 50, min_shared = 61)
  expect_true(all(keep_all))
  keep_few <- loci_ld_clump(ld_gt,
    thr_r2 = 0.2, size = 50, min_shared = 61, insufficient = "remove"
  )
  kept <- loci[keep_few, ]
  gaps <- tapply(kept$position, kept$chromosome, function(p) min(diff(p)))
  expect_true(all(gaps > 50000))
  expect_error(
    loci_ld_clump(ld_gt, thr_r2 = 0.2, min_shared = 1),
    "'min_shared' should be an integer of at least 2"
  )
  # imputed genotypes are used only if the gen_tibble is set to use them
  ld_gt <- gt_impute_simple(ld_gt, method = "mode")
  expect_identical(loci_ld_clump(ld_gt, thr_r2 = 0.2, size = 50), ref)
  gt_set_imputed(ld_gt, set = TRUE)
  imputed <- loci_ld_clump(ld_gt, thr_r2 = 0.2, size = 50)
  gt_set_imputed(ld_gt, set = FALSE)
  expect_false(identical(imputed, ref))
})

test_that("loci_ld_clump excludes loci and checks the cores", {
  # example-missing has missing genotypes at 367 of its 500 loci: keep the
  # 133 complete loci, without resaving the backing matrix
  bed_file <- system.file("extdata", "example-missing.bed", package = "bigsnpr")
  missing_gt <- gen_tibble(
    bed_file,
    backingfile = tempfile("missing_"),
    quiet = TRUE
  )
  complete_gt <- missing_gt %>%
    select_loci_if(loci_missingness(genotypes) == 0)
  keep_complete <- loci_ld_clump(complete_gt, thr_r2 = 0.2, size = 10)
  expect_length(keep_complete, 133)
  expect_error(
    loci_ld_clump(complete_gt, thr_r2 = 0.2, size = 10, n_cores = 0),
    "'ncores' should be at least 1"
  )
  # excluding the loci with missing genotypes gives the same clumps
  keep_excluded <- loci_ld_clump(missing_gt,
    thr_r2 = 0.2, size = 10,
    exclude = which(loci_missingness(missing_gt) > 0)
  )
  expect_identical(
    loci_names(missing_gt)[keep_excluded],
    loci_names(complete_gt)[keep_complete]
  )
})

test_that("ld clumping leaves no temporary FBMs in tempdir()", {
  # bigsnpr::snp_clumping() creates (and never removes) a temporary FBM for
  # each chromosome; their .bk files used to accumulate in tempdir(), where
  # a new tempfile() name could clash with them (tempfile() only checks the
  # name without the .bk)
  n_bk <- function() length(list.files(tempdir(), pattern = "\\.bk$"))
  bk_before <- n_bk()
  fbm_dir_before <- getOption("FBM.dir")
  keep <- loci_ld_clump(test_gt, thr_r2 = 0.2, return_id = TRUE)
  expect_identical(n_bk(), bk_before)
  # the option is restored
  expect_identical(getOption("FBM.dir"), fbm_dir_before)
})

test_that("with_private_fbm_dir() uses a fresh directory and deletes it", {
  fbm_dir_before <- getOption("FBM.dir")
  res <- with_private_fbm_dir({
    tmp_fbm <- bigstatsr::FBM(2, 2)
    list(
      dir = getOption("FBM.dir"),
      # (FBM() normalises the path of its backing file, with forward slashes
      # on Windows: both paths are normalised in the same way, while they
      # still exist)
      dir_norm = normalizePath(getOption("FBM.dir"), winslash = "/"),
      bk_dir_norm = normalizePath(dirname(tmp_fbm$backingfile),
        winslash = "/"
      ),
      bk = tmp_fbm$backingfile
    )
  })
  # the temporary FBM was created in the private directory
  expect_identical(res$bk_dir_norm, res$dir_norm)
  # and both are gone, and the option is restored
  expect_false(dir.exists(res$dir))
  expect_false(file.exists(res$bk))
  expect_identical(getOption("FBM.dir"), fbm_dir_before)
  # the option and the directory are cleaned up on error too
  expect_error(
    with_private_fbm_dir({
      err_dir <- getOption("FBM.dir")
      stop("boom")
    }),
    "boom"
  )
  expect_false(dir.exists(err_dir))
  expect_identical(getOption("FBM.dir"), fbm_dir_before)
  # the private directory is made inside the FBM.dir set by the user
  user_dir <- withr::local_tempdir()
  withr::local_options(FBM.dir = user_dir)
  private_dir <- with_private_fbm_dir(getOption("FBM.dir"))
  # (normalised, as dirname() changes the separators on Windows)
  expect_identical(
    normalizePath(dirname(private_dir), winslash = "/"),
    normalizePath(user_dir, winslash = "/")
  )
  expect_false(dir.exists(private_dir))
  expect_identical(getOption("FBM.dir"), user_dir)
  # a clear error if FBM.dir does not exist (and the option is unchanged)
  missing_dir <- file.path(user_dir, "not_a_dir")
  withr::local_options(FBM.dir = missing_dir)
  expect_error(
    with_private_fbm_dir(bigstatsr::FBM(2, 2)),
    "could not create a directory for temporary backing files"
  )
  expect_identical(getOption("FBM.dir"), missing_dir)
})

test_that("S statistic is the correct length", {
  s_stat <- loci_maf(test_gt)
  s_stat_subset <- s_stat[1:5]
  expect_error(
    loci_ld_clump(test_gt,
      thr_r2 = 0.2,
      return_id = TRUE,
      S = s_stat_subset
    ),
    "Length of 'S' must equal "
  )
  # but it works when using a full length S
  keep <- loci_ld_clump(test_gt, thr_r2 = 0.2, return_id = TRUE, S = s_stat)
  expect_true(all.equal(keep, c(1, 2, 3, 4, 6)) == TRUE)
})

test_that("loci_ld_clump returns the same as bigsnpr", {
  bedfile <- system.file("extdata/related/families.bed", package = "tidypopgen")
  rds <- bigsnpr::snp_readBed(bedfile, backingfile = tempfile())

  bigsnp <- bigsnpr::snp_attach(rds)

  gen_tbl <- gen_tibble(
    bedfile,
    quiet = TRUE,
    backingfile = tempfile(),
    valid_alleles = c("1", "2")
  )

  # also tests imputation
  bigsnp_imputed <- bigsnpr::snp_fastImputeSimple(
    bigsnp$genotypes,
    method = "mode"
  )
  bigsnp$genotypes <- bigsnp_imputed
  bigsnpr::snp_save(bigsnp)
  bigsnp_imputed <- bigsnpr::snp_attach(rds)

  gen_tbl_imputed <- gt_impute_simple(gen_tbl, method = "mode")
  gt_set_imputed(gen_tbl_imputed, TRUE)

  # set up bigsnpr
  G <- bigsnp_imputed$genotypes # nolint start
  POS <- bigsnp_imputed$map$physical.pos
  CHR <- bigsnp_imputed$map$chromosome

  ind.keep <- bigsnpr::snp_clumping( # nolint end
    G,
    infos.chr = CHR,
    infos.pos = POS,
    thr.r2 = 0.2
  )
  to_keep <- loci_ld_clump(gen_tbl_imputed, thr_r2 = 0.2)

  # convert our output to indices
  ind <- which(to_keep == TRUE)

  # check they match
  expect_true(all(ind.keep == ind))
})


test_that("loci_ld_clump error unsorted loci", {
  pop_b <- gen_tibble(
    system.file("extdata/pop_b.bed", package = "tidypopgen"),
    backingfile = tempfile(),
    quiet = TRUE
  )

  # now scramble the loci
  set.seed(123)
  random_order <- sample(1:17)
  show_loci(pop_b) <- pop_b %>%
    select_loci(all_of(random_order)) %>%
    show_loci()

  # impute (and use the imputed genotypes)
  pop_b_imputed <- gt_impute_simple(pop_b, method = "mode")
  gt_set_imputed(pop_b_imputed, set = TRUE)

  # ld
  expect_error(
    loci_ld_clump(pop_b_imputed, thr_r2 = 0.2),
    "Your loci have been resorted"
  )
  expect_false(identical(
    show_loci(pop_b_imputed),
    pop_b_imputed %>% show_loci() %>% arrange(chromosome, position)
  ))

  # reorder the loci
  show_loci(pop_b_imputed) <- pop_b_imputed %>%
    show_loci() %>%
    arrange(chromosome, position)

  # try again
  expect_equal(
    loci_ld_clump(pop_b_imputed, thr_r2 = 0.2),
    c(
      FALSE,
      TRUE,
      TRUE,
      FALSE,
      TRUE,
      TRUE,
      FALSE,
      FALSE,
      TRUE,
      TRUE,
      TRUE,
      TRUE,
      TRUE,
      TRUE,
      TRUE,
      TRUE,
      TRUE
    )
  )
  expect_true(identical(
    show_loci(pop_b_imputed),
    pop_b_imputed %>% show_loci() %>% arrange(chromosome, position)
  ))
})

test_that("loci order", {
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c"),
    population = c("pop1", "pop1", "pop2")
  )
  test_genotypes <- rbind(
    c(1, 2, 2, 0, 1, 2),
    c(0, 2, 2, 0, 0, 1),
    c(1, 2, 2, 0, 0, 1)
  )

  test_loci <- data.frame(
    name = paste0("rs", c(5, rep(1:4, 1), 6)),
    chromosome = c("chr2", "chr1", "chr1", "chr1", "chr1", "chr2"),
    position = c(23, 3, 5, 65, 343, 456),
    genetic_dist = as.double(rep(0, 6)),
    allele_ref = c("C", "A", "T", "C", "G", "T"),
    allele_alt = c("G", "T", "C", NA, "C", "A")
  )

  test_gt_new_order <- gen_tibble(
    x = test_genotypes,
    indiv_meta = test_indiv_meta,
    loci = test_loci,
    backingfile = tempfile(),
    quiet = TRUE
  )

  # clumping generates erro
  expect_error(
    loci_ld_clump(test_gt_new_order, thr_r2 = 0.2),
    "All SNPs in a chromosome should be adjacent in the loci table"
  )
  # reorder the loci
  show_loci(test_gt_new_order) <- test_gt_new_order %>%
    show_loci() %>%
    arrange(chromosome, position)

  # try again
  expect_error(
    loci_ld_clump(test_gt_new_order, thr_r2 = 0.2, return_id = TRUE),
    "Your loci have been resorted"
  )

  # TODO we need to resave the genotypes to the backing file in the right order

  # calculate the expected result
  #  keep <- loci_ld_clump(test_gt, thr_r2 = 0.2, return_id=TRUE) #nolint
  # compare
  #  expect_equal(keep, keep_reordered) #nolint
})

test_that("loci_ld_clump works on a grouped gt", {
  bedfile <- system.file("extdata/related/families.bed", package = "tidypopgen")
  rds <- bigsnpr::snp_readBed(bedfile, backingfile = tempfile())
  bigsnp <- bigsnpr::snp_attach(rds)

  gen_tbl <- gen_tibble(
    bedfile,
    quiet = TRUE,
    backingfile = tempfile(),
    valid_alleles = c("1", "2")
  )

  imputed_data <- gt_impute_simple(gen_tbl, method = "random")
  to_keep_ld_ungrouped <- loci_ld_clump(imputed_data, thr_r2 = 0.2, size = 10)

  gen_tbl$population <- rep(c("population_1", "population_2"), each = 6)
  gen_tbl <- gen_tbl %>% group_by(population)

  imputed_data <- gt_impute_simple(gen_tbl, method = "random")
  to_keep_ld_grouped <- loci_ld_clump(imputed_data, thr_r2 = 0.2, size = 10)

  # Removed loci are chosen at random, so we can't use expect equal
  # However, the same number of loci should be removed in both cases
  expect_equal(
    length(to_keep_ld_ungrouped == FALSE),
    length(to_keep_ld_grouped == FALSE)
  )
  expect_equal(
    length(to_keep_ld_ungrouped == TRUE),
    length(to_keep_ld_grouped == TRUE)
  )
})

test_that("n_cores can be set", {
  expect_true(getOption("bigstatsr.check.parallel.blas"))
  one_core <- loci_ld_clump(test_gt, thr_r2 = 0.2, n_cores = 1)
  two_core <- loci_ld_clump(test_gt, thr_r2 = 0.2, n_cores = 2)
  expect_equal(one_core, two_core)
  expect_true(getOption("bigstatsr.check.parallel.blas"))

  # test parallel blas is true on exit if function errors
  # scramble the loci
  set.seed(123)
  random_order <- sample(1:6)
  show_loci(test_gt) <- test_gt %>%
    select_loci(all_of(random_order)) %>%
    show_loci()
  # impute
  test_gt_imputed <- gt_impute_simple(test_gt, method = "mode")
  expect_error(
    loci_ld_clump(test_gt_imputed, thr_r2 = 0.2, n_cores = 2),
    "Your loci have been resorted"
  )
  expect_true(getOption("bigstatsr.check.parallel.blas"))
})
