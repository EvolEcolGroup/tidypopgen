skip_if_not_installed("admixtools")

test_genotypes <- rbind(
  c(1, 1, 0, 1, 1, 0),
  c(2, 1, 0, NA, 0, 0),
  c(2, NA, 0, 0, 1, 1),
  c(2, 0, 0, 2, 0, 0),
  c(1, 2, 0, 1, 2, 1),
  c(0, 0, 0, 0, NA, 2),
  c(0, 1, 2, 0, 1, NA)
)
test_indiv_meta <- data.frame(
  id = c("a", "b", "c", "d", "e", "f", "g"),
  population = c("pop1", "pop1", "pop2", "pop2", "pop1", "pop3", "pop3")
)
test_loci <- data.frame(
  name = paste0("rs", 1:6),
  chromosome = paste0("chr", c(1, 1, 1, 1, 2, 2)),
  position = as.integer(c(3, 5, 65, 343, 23, 456)),
  genetic_dist = as.double(rep(0, 6)),
  allele_ref = c("A", "T", "C", "G", "C", "T"),
  allele_alt = c("T", "C", NA, "C", "G", "A")
)

test_gt <- gen_tibble(
  x = test_genotypes,
  loci = test_loci,
  indiv_meta = test_indiv_meta,
  quiet = TRUE
)
test_gt <- test_gt %>% group_by(population)

test_that("extract f2 correctly", {
  # process the data with admixtools (note that we get some warnings)
  bed_file <- gt_as_plink(test_gt, file = tempfile("test_bed"))

  test_gt2 <- gen_tibble(bed_file, quiet = TRUE)

  # test af table
  # without adjusting pseudohaploids
  adm_aftable <- admixtools:::anygeno_to_aftable(
    bigsnpr::sub_bed(bed_file),
    adjust_pseudohaploid = FALSE,
    verbose = FALSE
  )
  gt_aftable <- gt_to_aftable(test_gt)

  expect_true(all.equal(
    adm_aftable$afs,
    gt_aftable$afs,
    check.attributes = FALSE
  ))
  expect_true(all.equal(
    adm_aftable$counts,
    gt_aftable$counts,
    check.attributes = FALSE
  ))
  expect_true(all.equal(
    adm_aftable$snpfile,
    gt_aftable$snpfile %>% select(-chr_int),
    check.attributes = FALSE
  ))

  # now adjusting the pseudohaploids
  adm_aftable <- admixtools:::anygeno_to_aftable(
    bigsnpr::sub_bed(bed_file),
    adjust_pseudohaploid = TRUE,
    verbose = FALSE
  )
  test_gt <- gt_pseudohaploid(test_gt)
  gt_aftable <- gt_to_aftable(test_gt)
  # expect_true(all.equal(adm_aftable, gt_aftable, check.attributes= FALSE)) #nolint

  expect_true(all.equal(
    adm_aftable$afs,
    gt_aftable$afs,
    check.attributes = FALSE
  ))
  expect_true(all.equal(
    adm_aftable$counts,
    gt_aftable$counts,
    check.attributes = FALSE
  ))
  expect_true(all.equal(
    adm_aftable$snpfile,
    gt_aftable$snpfile %>% select(-chr_int),
    check.attributes = FALSE
  ))

  adm_outdir <- file.path(tempdir(), "adm_f2")
  unlink(file.path(adm_outdir, "*"), recursive = TRUE)
  # we get a few warnings due to the small sample size
  suppressWarnings(admixtools::extract_f2(
    bigsnpr::sub_bed(bed_file),
    outdir = adm_outdir,
    verbose = FALSE
  ))
  expect_warning(
    adm_f2 <- admixtools::f2_from_precomp(adm_outdir, verbose = FALSE)
  ) # nolint
  # now try to do the same with gen_tibble
  gt_outdir <- file.path(tempdir(), "gt_f2")
  unlink(file.path(gt_outdir, "*"), recursive = TRUE)
  # we get same warning due to the small dataset
  # TODO can we capture the warnings above and then check that they are the same
  suppressWarnings(gt_extract_f2(test_gt, outdir = gt_outdir, quiet = TRUE))
  expect_warning(
    gt_f2 <- admixtools::f2_from_precomp(adm_outdir, verbose = FALSE)
  ) # nolint
  expect_true(all.equal(adm_f2, gt_f2))
})

test_that("n_cores can be set", {
  temp1 <- tempfile()
  temp2 <- tempfile()

  suppressWarnings(gt_extract_f2(test_gt,
    outdir = temp1,
    quiet = TRUE,
    n_cores = 1
  ))
  suppressWarnings(gt_extract_f2(test_gt,
    outdir = temp2,
    quiet = TRUE,
    n_cores = 2
  ))

  expect_equal(
    suppressWarnings(admixtools::f2_from_precomp(temp1, verbose = FALSE)),
    suppressWarnings(admixtools::f2_from_precomp(temp2, verbose = FALSE))
  )
})
