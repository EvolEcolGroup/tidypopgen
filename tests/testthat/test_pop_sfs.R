# Tests for pop_sfs() and pop_sfs_preview(), against the reference
# implementation of easySFS/dadi in helper_sfs_reference.R

# a whole population missing at a locus, and a locus with a single genotype
# in a population
sim <- ref_sfs_simulate(c(popA = 8, popB = 6, popC = 5),
  n_loci = 120,
  modify = function(genotypes, population) {
    genotypes[population == "popA", 10] <- NA
    genotypes[population == "popC", 20] <- c(1, NA, NA, NA, NA)
    genotypes
  }
)
pops <- c("popA", "popB", "popC")

test_that("pop_sfs matches easySFS (folded and unfolded)", {
  for (n_proj in list(c(10, 8, 6), c(16, 12, 10), c(2, 2, 2), c(1, 5, 3))) {
    for (folded in c(TRUE, FALSE)) {
      sfs <- pop_sfs(sim$gt,
        n_proj = n_proj, folded = folded,
        type = "array"
      )
      expect_named(sfs, pops)
      for (i in seq_along(pops)) {
        ref <- ref_sfs(sim$alt_copies, sim$ploidy, sim$population,
          pops = pops[i], proj = n_proj[i], polarized = !folded
        )
        expect_equal(sfs_values(sfs[[i]]), sfs_values(ref), tolerance = 1e-10)
        expect_identical(attr(sfs[[i]], "mask"), unname(attr(ref, "mask")))
        expect_identical(attr(sfs[[i]], "folded"), folded)
        expect_identical(
          attr(sfs[[i]], "n_proj"),
          stats::setNames(as.integer(n_proj[i]), pops[i])
        )
        # the spectrum sums to the number of loci used
        expect_equal(sum(sfs[[i]]), attr(sfs[[i]], "n_loci"))
        expect_equal(
          attr(sfs[[i]], "n_loci") + attr(sfs[[i]], "n_loci_dropped"),
          count_loci(sim$gt)
        )
      }
    }
  }
})

test_that("pop_sfs returns a tidy tibble", {
  sfs_tidy <- pop_sfs(sim$gt, n_proj = c(4, 6, 2))
  expect_s3_class(sfs_tidy, "tbl_df")
  expect_named(sfs_tidy, c("group", "allele_count", "n_sites", "masked"))
  expect_identical(sfs_tidy$group, rep(pops, c(5, 7, 3)))
  expect_identical(sfs_tidy$allele_count, c(0:4, 0:6, 0:2))
  expect_identical(attr(sfs_tidy, "n_proj"), c(popA = 4L, popB = 6L, popC = 2L))
  expect_true(attr(sfs_tidy, "folded"))
  sfs_array <- pop_sfs(sim$gt, n_proj = c(4, 6, 2), type = "array")
  expect_equal(
    sfs_tidy$n_sites,
    unlist(lapply(sfs_array, as.vector), use.names = FALSE)
  )
  expect_equal(
    sfs_tidy$masked,
    unlist(lapply(sfs_array, function(x) as.vector(attr(x, "mask"))),
      use.names = FALSE
    )
  )
  # tidy() of each sfs gives the same cells
  expect_equal(
    tidy(sfs_array$popB)$n_sites,
    sfs_tidy$n_sites[sfs_tidy$group == "popB"]
  )
  expect_named(tidy(sfs_array$popB), c("popB", "n_sites", "masked"))
  # a factor grouping variable keeps its levels
  gt_factor <- sim$gt %>%
    dplyr::ungroup() %>%
    dplyr::mutate(population = factor(population,
      levels = c("popC", "popA", "popB")
    )) %>%
    dplyr::group_by(population)
  sfs_factor <- pop_sfs(gt_factor, n_proj = c(popA = 4, popB = 6, popC = 2))
  expect_identical(levels(sfs_factor$group), c("popC", "popA", "popB"))
  expect_equal(
    sfs_factor$n_sites[sfs_factor$group == "popA"],
    sfs_tidy$n_sites[sfs_tidy$group == "popA"]
  )
})

test_that("pop_sfs works on ungrouped data", {
  all_pops <- rep("all", length(sim$population))
  ref <- ref_sfs(sim$alt_copies, sim$ploidy, all_pops,
    pops = "all", proj = 20
  )
  sfs <- pop_sfs(dplyr::ungroup(sim$gt), n_proj = 20, type = "array")
  expect_s3_class(sfs, "sfs")
  expect_equal(sfs_values(sfs), sfs_values(ref), tolerance = 1e-10)
  sfs_vctrs <- pop_sfs(sim$gt$genotypes, n_proj = 20, type = "array")
  expect_identical(sfs_vctrs, sfs)
  sfs_tidy <- pop_sfs(dplyr::ungroup(sim$gt), n_proj = 20)
  expect_named(sfs_tidy, c("allele_count", "n_sites", "masked"))
  expect_equal(sfs_tidy$n_sites, as.vector(sfs))
})

test_that("pop_sfs without missing data is the count of alternate copies", {
  sim_full <- ref_sfs_simulate(c(popA = 5, popB = 4),
    n_loci = 50,
    max_missing = 0, seed = 2
  )
  sfs <- pop_sfs(sim_full$gt,
    n_proj = c(10, 8), folded = FALSE,
    type = "array"
  )
  alt_a <- colSums(sim_full$alt_copies[sim_full$population == "popA", ])
  expect_equal(as.vector(sfs$popA), tabulate(alt_a + 1, nbins = 11))
  # the folded spectrum counts the minor allele
  sfs_folded <- pop_sfs(sim_full$gt, n_proj = c(10, 8), type = "array")
  minor_a <- pmin(alt_a, 10 - alt_a)
  expected <- tabulate(minor_a + 1, nbins = 11)
  # the cell with half the copies is averaged with itself, i.e. unchanged
  expect_equal(as.vector(sfs_folded$popA), expected)
  expect_identical(
    as.vector(attr(sfs_folded$popA, "mask")),
    c(TRUE, rep(FALSE, 5), rep(TRUE, 5))
  )
})

test_that("pop_sfs does not depend on the number of cores or block size", {
  sfs_1 <- pop_sfs(sim$gt, n_proj = c(10, 8, 6), n_cores = 1, type = "array")
  sfs_2 <- pop_sfs(sim$gt,
    n_proj = c(10, 8, 6), n_cores = 2, block_size = 7,
    type = "array"
  )
  expect_identical(sfs_1, sfs_2)
})

test_that("pop_sfs checks n_proj", {
  expect_error(pop_sfs(sim$gt, n_proj = c(4, 4)), "one value per group")
  expect_error(pop_sfs(sim$gt, n_proj = 17), "popA \\(16\\)")
  expect_error(pop_sfs(sim$gt, n_proj = 0), "at least 1")
  expect_error(pop_sfs(sim$gt, n_proj = 2.5), "whole numbers")
  expect_error(
    pop_sfs(sim$gt, n_proj = c(popA = 2, popB = 2, popD = 2)),
    "not groups: popD"
  )
  expect_error(
    pop_sfs(sim$gt, n_proj = c(popA = 2, popB = 2)),
    "no value for: popC"
  )
  expect_error(
    pop_sfs(sim$gt, n_proj = c(popA = 2, popB = 2, popC = 2, popA = 4)),
    "duplicate names: popA"
  )
  expect_error(pop_sfs(sim$gt, n_proj = 2, folded = NA), "TRUE or FALSE")
  # named values can be in any order
  expect_identical(
    pop_sfs(sim$gt, n_proj = c(popC = 2, popA = 4, popB = 6)),
    pop_sfs(sim$gt, n_proj = c(4, 6, 2))
  )
})

test_that("pop_sfs_preview matches easySFS --preview", {
  preview <- pop_sfs_preview(sim$gt)
  expect_named(
    preview,
    c("group", "n_proj", "n_seg_sites", "n_loci", "best")
  )
  max_n <- c(popA = 16, popB = 12, popC = 10)
  for (pop in pops) {
    this_preview <- preview[preview$group == pop, ]
    expect_identical(this_preview$n_proj, 2:max_n[[pop]])
    ref_seg <- vapply(2:max_n[[pop]], function(p) {
      ref_sfs_seg_sites(ref_sfs(sim$alt_copies, sim$ploidy, sim$population,
        pops = pop, proj = p
      ))
    }, numeric(1))
    expect_equal(this_preview$n_seg_sites, ref_seg, tolerance = 1e-10)
    # the same for the unfolded spectrum
    ref_seg_unfolded <- vapply(2:max_n[[pop]], function(p) {
      ref_sfs_seg_sites(ref_sfs(sim$alt_copies, sim$ploidy, sim$population,
        pops = pop, proj = p, polarized = TRUE
      ))
    }, numeric(1))
    expect_equal(this_preview$n_seg_sites, ref_seg_unfolded, tolerance = 1e-10)
    # number of loci with enough called copies
    called <- ref_sfs_calls(sim$alt_copies, sim$ploidy, sim$population, pop)
    expect_equal(
      this_preview$n_loci,
      vapply(2:max_n[[pop]], function(p) {
        sum(called[[1]][, "called"] >= p)
      }, numeric(1))
    )
    # the best projection maximises the number of segregating sites
    expect_identical(sum(this_preview$best), 1L)
    is_max <- this_preview$n_seg_sites == max(this_preview$n_seg_sites)
    expect_identical(
      this_preview$n_proj[this_preview$best],
      max(this_preview$n_proj[is_max])
    )
  }
  # ungrouped data
  preview_all <- pop_sfs_preview(dplyr::ungroup(sim$gt))
  expect_named(preview_all, c("n_proj", "n_seg_sites", "n_loci", "best"))
  expect_identical(preview_all$n_proj, 2:38)
})

test_that("pop_sfs chooses the projection with pop_sfs_preview by default", {
  best <- pop_sfs_preview(sim$gt) %>% dplyr::filter(best)
  expect_message(
    sfs_default <- pop_sfs(sim$gt),
    "maximises its number of segregating sites"
  )
  expect_identical(
    attr(sfs_default, "n_proj"),
    stats::setNames(best$n_proj, best$group)
  )
  expect_identical(
    sfs_default,
    pop_sfs(sim$gt, n_proj = stats::setNames(best$n_proj, best$group))
  )
})

test_that("pop_sfs works with polyploid and mixed-ploidy data", {
  sim_poly <- ref_sfs_simulate(c(popA = 4, popB = 5),
    n_loci = 40,
    ploidy = c(4, 4, 2, 4, 3, 4, 4, 2, 4), seed = 3
  )
  sfs <- pop_sfs(sim_poly$gt, n_proj = c(8, 10), type = "array")
  for (i in 1:2) {
    ref <- ref_sfs(sim_poly$alt_copies, sim_poly$ploidy, sim_poly$population,
      pops = c("popA", "popB")[i], proj = c(8, 10)[i]
    )
    expect_equal(sfs_values(sfs[[i]]), sfs_values(ref), tolerance = 1e-10)
  }
  expect_error(pop_sfs(sim_poly$gt, n_proj = c(15, 2)), "popA \\(14\\)")
})

test_that("pop_sfs counts one allele copy for pseudohaploids", {
  test_genotypes <- rbind(
    c(0, 2, NA, 2, 0),
    c(1, 2, 0, 1, 1),
    c(0, NA, 2, 2, 0),
    c(2, 0, 0, NA, 2),
    c(1, 1, 2, 0, 0)
  )
  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = data.frame(
      name = paste0("rs", 1:5), chromosome = 1, position = 1:5,
      genetic_dist = 0, allele_ref = "A", allele_alt = "T"
    ),
    indiv_meta = data.frame(
      id = letters[1:5],
      population = c("popA", "popA", "popA", "popB", "popB")
    ),
    quiet = TRUE
  ) %>%
    gt_pseudohaploid(test_n_loci = NULL) %>%
    dplyr::group_by(population)
  ploidy <- indiv_ploidy(test_gt)
  expect_identical(ploidy, c(1, 2, 1, 1, 2))
  # pseudohaploids carry a single copy (their 0/2 dosage halved)
  alt_copies <- test_genotypes / (3 - ploidy)
  sfs <- pop_sfs(test_gt, n_proj = c(3, 2), folded = FALSE, type = "array")
  for (i in 1:2) {
    ref <- ref_sfs(alt_copies, ploidy, test_gt$population,
      pops = c("popA", "popB")[i], proj = c(3, 2)[i], polarized = TRUE
    )
    expect_equal(sfs_values(sfs[[i]]), sfs_values(ref), tolerance = 1e-10)
  }
  expect_error(pop_sfs(test_gt, n_proj = c(5, 2)), "popA \\(4\\)")
})

test_that("pop_sfs needs integer genotypes", {
  test_gt <- load_example_gt("grouped_gen_tbl")
  attr(test_gt$genotypes, "fbm")$code256 <- bigsnpr::CODE_DOSAGE
  expect_error(pop_sfs(test_gt, n_proj = 2), "cannot be computed from dosages")
  attr(test_gt$genotypes, "fbm")$code256 <- bigsnpr::CODE_012
})

test_that("pop_sfs drops empty groups", {
  gt_empty <- sim$gt %>%
    dplyr::ungroup() %>%
    dplyr::mutate(population = factor(population,
      levels = c("popA", "popB", "popC", "popD")
    )) %>%
    dplyr::group_by(population, .drop = FALSE)
  expect_message(
    sfs <- pop_sfs(gt_empty, n_proj = 2),
    "empty groups are dropped: popD"
  )
  expect_identical(
    as.character(unique(sfs$group)),
    c("popA", "popB", "popC")
  )
})

test_that("pop_sfs uses one column per grouping variable", {
  gt_two <- sim$gt %>%
    dplyr::ungroup() %>%
    dplyr::mutate(sex = rep(c("F", "M"), length.out = dplyr::n())) %>%
    dplyr::group_by(population, sex)
  sfs <- pop_sfs(gt_two, n_proj = 2)
  expect_named(sfs, c("population", "sex", "allele_count", "n_sites", "masked"))
  sfs_array <- pop_sfs(gt_two, n_proj = 2, type = "array")
  expect_named(sfs_array, c(
    "popA.F", "popA.M", "popB.F", "popB.M", "popC.F",
    "popC.M"
  ))
  # the labels of the groups (pasted with ".") must be distinct
  gt_same_label <- sim$gt %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      v1 = dplyr::if_else(population == "popA", "a.b", "a"),
      v2 = dplyr::if_else(population == "popA", "c", "b.c")
    ) %>%
    dplyr::group_by(v1, v2)
  expect_error(pop_sfs(gt_same_label, n_proj = 2), "same label.*: a.b.c;")
  # a grouping variable cannot have the name of another column
  gt_clash <- gt_two %>% dplyr::rename(n_sites = "sex")
  expect_error(pop_sfs(gt_clash, n_proj = 2), "other columns of the output")
  expect_error(
    bootstrap_pop_sfs(gt_two %>% dplyr::rename(replicate = "sex"),
      n_bootstrap = 2, block_size = 20, n_proj = 2
    ),
    "\\(replicate\\)"
  )
})
