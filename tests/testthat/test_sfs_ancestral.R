# Tests for the polarisation of unfolded spectra with the ancestral allele
# (the ancestral argument of pop_sfs(), pairwise_pop_sfs(), nwise_pop_sfs(),
# pop_sfs_preview() and write_sfs())

sim <- ref_sfs_simulate(c(popA = 6, popB = 5, popC = 4),
  n_loci = 80,
  seed = 41
)
pops <- c("popA", "popB", "popC")
loci <- show_loci(sim$gt)
# the reference allele is A and the alternate T at all loci; the ancestral
# allele is either of them (in either case), unknown, or another allele
withr::with_seed(42, {
  ancestral <- sample(c("A", "T", "a", "t", NA, "G"), nrow(loci),
    replace = TRUE, prob = c(0.35, 0.35, 0.1, 0.1, 0.05, 0.05)
  )
})
n_unknown <- sum(is.na(ancestral) | ancestral == "G")
polarised <- ref_sfs_polarise(sim$alt_copies, sim$ploidy, ancestral,
  allele_ref = loci$allele_ref, allele_alt = loci$allele_alt
)

test_that("unfolded spectra are polarised with the ancestral allele", {
  for (n_proj in list(c(8, 6, 4), c(12, 10, 8), c(3, 4, 2))) {
    sfs_1d <- pop_sfs(sim$gt,
      n_proj = n_proj, folded = FALSE,
      ancestral = ancestral, type = "array"
    )
    sfs_2d <- pairwise_pop_sfs(sim$gt,
      n_proj = n_proj, folded = FALSE,
      ancestral = ancestral, type = "array"
    )
    sfs_3d <- nwise_pop_sfs(sim$gt,
      n_proj = n_proj, folded = FALSE,
      ancestral = ancestral, type = "array"
    )
    for (i in seq_along(pops)) {
      ref <- ref_sfs(polarised, sim$ploidy, sim$population,
        pops = pops[i], proj = n_proj[i], polarized = TRUE
      )
      expect_equal(sfs_values(sfs_1d[[i]]), sfs_values(ref), tolerance = 1e-10)
    }
    pairs <- utils::combn(3, 2)
    for (i in seq_len(ncol(pairs))) {
      ref <- ref_sfs(polarised, sim$ploidy, sim$population,
        pops = pops[pairs[, i]], proj = n_proj[pairs[, i]], polarized = TRUE
      )
      expect_equal(sfs_values(sfs_2d[[i]]), sfs_values(ref), tolerance = 1e-10)
    }
    ref <- ref_sfs(polarised, sim$ploidy, sim$population,
      pops = pops, proj = n_proj, polarized = TRUE
    )
    expect_equal(sfs_values(sfs_3d), sfs_values(ref), tolerance = 1e-10)
    # loci are either used, dropped for missing data, or unpolarised
    expect_identical(attr(sfs_3d, "n_loci_unpolarised"), as.numeric(n_unknown))
    expect_equal(
      attr(sfs_3d, "n_loci") + attr(sfs_3d, "n_loci_dropped") +
        attr(sfs_3d, "n_loci_unpolarised"),
      count_loci(sim$gt)
    )
  }
})

test_that("ancestral can be a column of the loci table", {
  gt_anc <- sim$gt
  show_loci(gt_anc) <- show_loci(gt_anc) %>%
    dplyr::mutate(ancestral_allele = ancestral)
  expect_identical(
    pop_sfs(gt_anc,
      n_proj = 6, folded = FALSE,
      ancestral = "ancestral_allele"
    ),
    pop_sfs(sim$gt, n_proj = 6, folded = FALSE, ancestral = ancestral)
  )
  expect_identical(
    nwise_pop_sfs(gt_anc,
      n_proj = 3, folded = FALSE,
      ancestral = "ancestral_allele", type = "array"
    ),
    nwise_pop_sfs(sim$gt,
      n_proj = 3, folded = FALSE, ancestral = ancestral,
      type = "array"
    )
  )
})

test_that("the reference or alternate allele as ancestral", {
  sfs_null <- nwise_pop_sfs(sim$gt,
    n_proj = c(6, 4, 3), folded = FALSE,
    type = "array"
  )
  # the reference allele as ancestral is the default
  expect_identical(
    nwise_pop_sfs(sim$gt,
      n_proj = c(6, 4, 3), folded = FALSE,
      ancestral = loci$allele_ref, type = "array"
    ),
    sfs_null
  )
  # the alternate allele as ancestral reverses the spectrum along all axes
  sfs_alt <- nwise_pop_sfs(sim$gt,
    n_proj = c(6, 4, 3), folded = FALSE,
    ancestral = loci$allele_alt, type = "array"
  )
  expect_equal(as.vector(sfs_alt), rev(as.vector(sfs_null)))
  expect_identical(attr(sfs_alt, "mask"), attr(sfs_null, "mask"))
})

test_that("ancestral is ignored for folded spectra", {
  expect_warning(
    sfs_folded <- pairwise_pop_sfs(sim$gt,
      n_proj = 4,
      ancestral = ancestral
    ),
    "ignored for folded spectra"
  )
  expect_identical(sfs_folded, pairwise_pop_sfs(sim$gt, n_proj = 4))
})

test_that("ancestral is checked", {
  expect_error(
    pop_sfs(sim$gt, n_proj = 4, folded = FALSE, ancestral = c("A", "T")),
    "one value per locus \\(80\\)"
  )
  expect_error(
    pop_sfs(sim$gt, n_proj = 4, folded = FALSE, ancestral = "not_a_column"),
    "one value per locus"
  )
  expect_error(
    pop_sfs(sim$gt,
      n_proj = 4, folded = FALSE,
      ancestral = seq_len(nrow(loci))
    ),
    "character vector"
  )
})

test_that("pop_sfs_preview leaves out the loci that cannot be polarised", {
  preview <- pop_sfs_preview(sim$gt, ancestral = ancestral)
  for (pop in pops) {
    this_preview <- preview[preview$group == pop, ]
    ref_seg <- vapply(this_preview$n_proj, function(p) {
      ref_sfs_seg_sites(ref_sfs(polarised, sim$ploidy, sim$population,
        pops = pop, proj = p, polarized = TRUE
      ))
    }, numeric(1))
    expect_equal(this_preview$n_seg_sites, ref_seg, tolerance = 1e-10)
    called <- ref_sfs_calls(polarised, sim$ploidy, sim$population, pop)
    expect_equal(this_preview$n_loci[1], sum(called[[1]][, "called"] >= 2))
  }
  # the default projections use this preview
  best <- preview %>% dplyr::filter(best)
  expect_message(
    sfs_default <- pop_sfs(sim$gt, folded = FALSE, ancestral = ancestral),
    "maximises"
  )
  expect_identical(
    attr(sfs_default, "n_proj"),
    stats::setNames(best$n_proj, best$group)
  )
})

test_that("polarised spectra do not depend on the number of cores", {
  expect_identical(
    pairwise_pop_sfs(sim$gt,
      n_proj = 6, folded = FALSE, ancestral = ancestral,
      n_cores = 1, type = "array"
    ),
    pairwise_pop_sfs(sim$gt,
      n_proj = 6, folded = FALSE, ancestral = ancestral,
      n_cores = 2, block_size = 7, type = "array"
    )
  )
})

test_that("monomorphic loci are polarised with their reference allele", {
  test_gt <- load_example_gt("grouped_gen_tbl")
  # rs3 is monomorphic (no alternate allele): its reference allele (C) can be
  # ancestral, but not another allele
  expect_true(is.na(show_loci(test_gt)$allele_alt[3]))
  sfs_ref <- pop_sfs(test_gt,
    n_proj = 2, folded = FALSE, type = "array",
    ancestral = c("A", "T", "C", "G", "C", "T")
  )
  expect_identical(attr(sfs_ref$pop1, "n_loci_unpolarised"), 0)
  sfs_other <- pop_sfs(test_gt,
    n_proj = 2, folded = FALSE, type = "array",
    ancestral = c("A", "T", "G", "G", "C", "T")
  )
  expect_identical(attr(sfs_other$pop1, "n_loci_unpolarised"), 1)
  expect_output(print(sfs_other$pop1), "without ancestral allele: 1")
})

test_that("write_sfs polarises the unfolded spectra of a gen_tibble", {
  fsc_files <- write_sfs(sim$gt,
    format = "fastsimcoal2",
    out_dir = withr::local_tempdir(), n_proj = c(6, 4, 3),
    prefix = "test", folded = FALSE, ancestral = ancestral
  )
  expect_true("test_jointDAFpop1_0.obs" %in% basename(fsc_files))
  dadi_dir <- withr::local_tempdir()
  write_sfs(sim$gt,
    out_dir = dadi_dir, n_proj = c(6, 4, 3),
    folded = FALSE, ancestral = ancestral
  )
  expected <- withr::local_tempfile()
  write_sfs(
    nwise_pop_sfs(sim$gt,
      n_proj = c(6, 4, 3), folded = FALSE,
      ancestral = ancestral, type = "array"
    ),
    expected
  )
  expect_identical(
    readLines(file.path(dadi_dir, "popA-popB-popC.sfs")),
    readLines(expected)
  )
})
