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

  # impute
  pop_b_imputed <- gt_impute_simple(pop_b, method = "mode")

  # ld
  expect_error(
    loci_ld_clump(pop_b_imputed, thr_r2 = 0.2),
    "Your loci have been resorted"
  )
  expect_false(identical(
    show_loci(pop_b_imputed),
    pop_b_imputed %>% show_loci() %>% arrange(chr_int, position)
  ))

  # reorder the loci
  show_loci(pop_b_imputed) <- pop_b_imputed %>%
    show_loci() %>%
    arrange(chr_int, position)

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
    pop_b_imputed %>% show_loci() %>% arrange(chr_int, position)
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
    arrange(chr_int, position)

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
