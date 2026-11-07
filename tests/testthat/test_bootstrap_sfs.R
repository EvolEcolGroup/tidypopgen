# Tests for the block bootstrap of site frequency spectra
# (bootstrap_pop_sfs(), bootstrap_pairwise_pop_sfs(),
# bootstrap_nwise_pop_sfs() and the n_bootstrap argument of write_sfs())

sim <- ref_sfs_simulate(c(popA = 6, popB = 5, popC = 4),
  n_loci = 120,
  seed = 61
)
pops <- c("popA", "popB", "popC")
n_proj <- c(popA = 8L, popB = 6L, popC = 4L)

# a grouped gen_tibble with the loci of the blocks drawn in a replicate (each
# block repeated as many times as it was drawn)
resampled_gt <- function(sim, block_id, weights) {
  loci <- unlist(lapply(seq_along(weights), function(b) {
    rep(which(block_id == b), weights[b])
  }))
  gen_tibble(
    x = sim$alt_copies[, loci, drop = FALSE],
    loci = data.frame(
      name = paste0("rs", seq_along(loci)), chromosome = 1,
      position = seq_along(loci), genetic_dist = 0,
      allele_ref = "A", allele_alt = "T"
    ),
    indiv_meta = data.frame(
      id = paste0("ind", seq_along(sim$population)),
      population = sim$population
    ),
    backingfile = tempfile(),
    quiet = TRUE
  ) %>% dplyr::group_by(.data$population)
}

test_that("each replicate is the spectrum of the resampled blocks", {
  block_id <- .jackknife_block_ids(show_loci(sim$gt),
    block_size = 10,
    size_unit = "snp"
  )
  for (folded in c(TRUE, FALSE)) {
    withr::with_seed(1, {
      boot_1d <- bootstrap_pop_sfs(sim$gt,
        n_bootstrap = 3, block_size = 10,
        n_proj = n_proj, folded = folded, type = "array"
      )
    })
    withr::with_seed(2, {
      boot_2d <- bootstrap_pairwise_pop_sfs(sim$gt,
        n_bootstrap = 3,
        block_size = 10, n_proj = n_proj, folded = folded, type = "array"
      )
    })
    withr::with_seed(3, {
      boot_3d <- bootstrap_nwise_pop_sfs(sim$gt,
        n_bootstrap = 3, block_size = 10,
        n_proj = n_proj, folded = folded, type = "array"
      )
    })
    for (i_rep in 1:3) {
      gt_1d <- resampled_gt(
        sim, block_id,
        attr(boot_1d, "bootstrap_blocks")[i_rep, ]
      )
      expect_equal(
        boot_1d[[i_rep]],
        pop_sfs(gt_1d, n_proj = n_proj, folded = folded, type = "array"),
        tolerance = 1e-10
      )
      gt_2d <- resampled_gt(
        sim, block_id,
        attr(boot_2d, "bootstrap_blocks")[i_rep, ]
      )
      expect_equal(
        boot_2d[[i_rep]],
        pairwise_pop_sfs(gt_2d,
          n_proj = n_proj, folded = folded,
          type = "array"
        ),
        tolerance = 1e-10
      )
      gt_3d <- resampled_gt(
        sim, block_id,
        attr(boot_3d, "bootstrap_blocks")[i_rep, ]
      )
      expect_equal(
        boot_3d[[i_rep]],
        nwise_pop_sfs(gt_3d, n_proj = n_proj, folded = folded, type = "array"),
        tolerance = 1e-10
      )
    }
  }
})

test_that("the spectra of the blocks add up to the spectrum of all loci", {
  data <- .sfs_data(sim$gt)
  block_id <- .jackknife_block_ids(show_loci(sim$gt),
    block_size = 7,
    size_unit = "snp"
  )
  spectra <- list(1, c(1, 3), 1:3)
  configs <- .sfs_config_counts(data,
    spectra = spectra,
    min_n = n_proj, n_cores = 1, block_size = 50, block_ids = block_id
  )
  # a single replicate in which every block is drawn once
  all_blocks <- .sfs_project(configs, spectra, n_proj,
    folded = TRUE,
    n_cores = 1, weights = matrix(1, nrow = 1, ncol = max(block_id))
  )
  expect_equal(
    lapply(all_blocks, .sfs_replicate),
    list(
      pop_sfs(sim$gt, n_proj = n_proj, type = "array")$popA,
      pairwise_pop_sfs(sim$gt, n_proj = n_proj, type = "array")$popA.popC,
      nwise_pop_sfs(sim$gt, n_proj = n_proj, type = "array")
    ),
    tolerance = 1e-10
  )
})

test_that("bootstrap replicates are reproducible", {
  boot_1 <- withr::with_seed(10, {
    bootstrap_pairwise_pop_sfs(sim$gt,
      n_bootstrap = 5, block_size = 10,
      n_proj = n_proj, n_cores = 1
    )
  })
  boot_2 <- withr::with_seed(10, {
    bootstrap_pairwise_pop_sfs(sim$gt,
      n_bootstrap = 5, block_size = 10,
      n_proj = n_proj, n_cores = 2
    )
  })
  expect_identical(boot_1, boot_2)
  # each replicate draws as many blocks as there are
  blocks <- attr(boot_1, "bootstrap_blocks")
  expect_identical(dim(blocks), c(5L, 12L))
  expect_true(all(rowSums(blocks) == 12))
})

test_that("bootstrap replicates are returned as tidy tibbles", {
  boot <- withr::with_seed(11, {
    bootstrap_pop_sfs(sim$gt, n_bootstrap = 4, block_size = 20, n_proj = n_proj)
  })
  expect_named(boot, c(
    "group", "replicate", "allele_count", "n_sites",
    "masked"
  ))
  expect_identical(nrow(boot), 4L * sum(n_proj + 1L))
  boot_array <- withr::with_seed(11, {
    bootstrap_pop_sfs(sim$gt,
      n_bootstrap = 4, block_size = 20, n_proj = n_proj,
      type = "array"
    )
  })
  expect_identical(
    attr(boot, "bootstrap_blocks"),
    attr(boot_array, "bootstrap_blocks")
  )
  expect_identical(attr(boot, "n_proj"), n_proj)
  # the cells of replicate 3 of popB
  rep_3 <- boot[boot$group == "popB" & boot$replicate == 3, ]
  expect_equal(rep_3$n_sites, as.vector(boot_array[[3]]$popB))
  expect_identical(rep_3$masked, as.vector(attr(boot_array[[3]]$popB, "mask")))

  # pairwise, with the second population changing fastest
  boot_2d <- withr::with_seed(12, {
    bootstrap_pairwise_pop_sfs(sim$gt,
      n_bootstrap = 2, block_size = 20, n_proj = n_proj
    )
  })
  expect_named(boot_2d, c(
    "population_1", "population_2", "replicate",
    "allele_count_1", "allele_count_2", "n_sites", "masked"
  ))
  boot_2d_array <- withr::with_seed(12, {
    bootstrap_pairwise_pop_sfs(sim$gt,
      n_bootstrap = 2, block_size = 20,
      n_proj = n_proj, type = "array"
    )
  })
  is_ac_2 <- boot_2d$population_1 == "popA" &
    boot_2d$population_2 == "popC" & boot_2d$replicate == 2
  ac <- boot_2d[is_ac_2, ]
  expect_equal(ac$n_sites, as.vector(t(boot_2d_array[[2]]$popA.popC)))

  # nwise, with the first population changing fastest
  boot_3d <- withr::with_seed(13, {
    bootstrap_nwise_pop_sfs(sim$gt,
      n_bootstrap = 2, block_size = 20, n_proj = c(4, 4, 2),
      fastest = "first"
    )
  })
  expect_named(boot_3d, c("replicate", pops, "n_sites", "masked"))
  expect_identical(boot_3d$popA[1:5], c(0:4))
  # the cells of each replicate are those of its array (in column-major
  # order, as the first population changes fastest)
  boot_3d_array <- withr::with_seed(13, {
    bootstrap_nwise_pop_sfs(sim$gt,
      n_bootstrap = 2, block_size = 20, n_proj = c(4, 4, 2),
      type = "array"
    )
  })
  for (i_rep in 1:2) {
    expect_equal(
      boot_3d$n_sites[boot_3d$replicate == i_rep],
      as.vector(boot_3d_array[[i_rep]])
    )
  }

  # ungrouped data
  boot_all <- withr::with_seed(14, {
    bootstrap_pop_sfs(dplyr::ungroup(sim$gt),
      n_bootstrap = 2, block_size = 20, n_proj = 10
    )
  })
  expect_named(boot_all, c("replicate", "allele_count", "n_sites", "masked"))
})

test_that("blocks can be defined in base pairs or as single SNPs", {
  # loci are at positions 1 to 120, so blocks of 25 bp give 5 blocks
  boot_bp <- withr::with_seed(15, {
    bootstrap_pop_sfs(sim$gt,
      n_bootstrap = 2, block_size = 25,
      size_unit = "bp", n_proj = n_proj
    )
  })
  expect_identical(ncol(attr(boot_bp, "bootstrap_blocks")), 5L)
  boot_snp <- withr::with_seed(16, {
    bootstrap_pop_sfs(sim$gt, n_bootstrap = 2, block_size = 1, n_proj = n_proj)
  })
  expect_identical(ncol(attr(boot_snp, "bootstrap_blocks")), 120L)
})

test_that("projections are chosen with all loci, and ancestral is used", {
  best <- pop_sfs_preview(sim$gt) %>% dplyr::filter(best)
  expect_message(
    boot <- withr::with_seed(17, {
      bootstrap_pop_sfs(sim$gt, n_bootstrap = 2, block_size = 20)
    }),
    "maximises"
  )
  expect_identical(
    attr(boot, "n_proj"),
    stats::setNames(best$n_proj, best$group)
  )
  # all loci polarised with the alternate allele reverse each replicate
  boot_ref <- withr::with_seed(18, {
    bootstrap_nwise_pop_sfs(sim$gt,
      n_bootstrap = 2, block_size = 20,
      n_proj = n_proj, folded = FALSE, type = "array"
    )
  })
  boot_alt <- withr::with_seed(18, {
    bootstrap_nwise_pop_sfs(sim$gt,
      n_bootstrap = 2, block_size = 20,
      n_proj = n_proj, folded = FALSE, type = "array",
      ancestral = show_loci(sim$gt)$allele_alt
    )
  })
  expect_equal(as.vector(boot_alt[[2]]), rev(as.vector(boot_ref[[2]])))
})

test_that("bootstrap arguments are checked", {
  expect_error(
    bootstrap_pop_sfs(sim$gt, n_proj = 2),
    "`block_size` is needed"
  )
  expect_error(
    bootstrap_pop_sfs(sim$gt, n_bootstrap = 0, block_size = 10, n_proj = 2),
    "positive whole number"
  )
  expect_error(
    bootstrap_pop_sfs(sim$gt, block_size = 1000, n_proj = 2),
    "at least 2 blocks"
  )
  expect_error(
    bootstrap_pairwise_pop_sfs(dplyr::ungroup(sim$gt),
      block_size = 10,
      n_proj = 2
    ),
    "grouped gen_tibble"
  )
  expect_error(
    bootstrap_nwise_pop_sfs(sim$gt,
      n_bootstrap = 100, block_size = 10,
      n_proj = n_proj, max_cells = 1000
    ),
    "more than `max_cells`"
  )
  # the number of block weights (replicates times blocks) is capped, and
  # write_sfs() stops before writing any file
  withr::local_options(tidypopgen.sfs_max_bootstrap_weights = 100)
  expect_error(
    bootstrap_pop_sfs(sim$gt,
      n_bootstrap = 10, block_size = 1, n_proj = 2
    ),
    "1200 weights \\(10 replicates of 120 blocks\\)"
  )
  out_dir <- withr::local_tempdir()
  expect_error(
    write_sfs(sim$gt,
      out_dir = out_dir, n_proj = 2,
      n_bootstrap = 10, bootstrap_block_size = 1
    ),
    "too much memory"
  )
  expect_length(list.files(out_dir), 0)
})

test_that("write_sfs writes the bootstrap replicates of a gen_tibble", {
  out_dir <- withr::local_tempdir()
  files <- withr::with_seed(20, {
    write_sfs(sim$gt,
      format = "fastsimcoal2",
      out_dir = out_dir, n_proj = n_proj, prefix = "test",
      n_bootstrap = 3, bootstrap_block_size = 20
    )
  })
  expect_identical(length(files), 4L * 7L)
  expect_true(all(file.exists(files)))
  expect_true(all(dir.exists(file.path(out_dir, paste0("bootstrap_", 1:3)))))
  # the files of all loci are unchanged
  out_dir_no_boot <- withr::local_tempdir()
  write_sfs(sim$gt,
    format = "fastsimcoal2",
    out_dir = out_dir_no_boot, n_proj = n_proj,
    prefix = "test"
  )
  for (f in list.files(out_dir_no_boot, recursive = TRUE)) {
    expect_identical(
      readLines(file.path(out_dir, f)),
      readLines(file.path(out_dir_no_boot, f))
    )
  }
  # the replicates are those of bootstrap_nwise_pop_sfs() with the same seed,
  # in both formats
  boot_3d <- withr::with_seed(20, {
    bootstrap_nwise_pop_sfs(sim$gt,
      n_bootstrap = 3, block_size = 20,
      n_proj = n_proj, type = "array"
    )
  })
  dadi_dir <- withr::local_tempdir()
  withr::with_seed(20, {
    write_sfs(sim$gt,
      out_dir = dadi_dir, n_proj = n_proj,
      n_bootstrap = 3, bootstrap_block_size = 20
    )
  })
  for (i_rep in 1:3) {
    expected <- withr::local_tempfile()
    write_sfs(boot_3d[[i_rep]], expected,
      format = "fastsimcoal2_multidimensional"
    )
    expect_identical(
      readLines(file.path(
        out_dir, paste0("bootstrap_", i_rep), "test_MSFS.obs"
      )),
      readLines(expected)
    )
    write_sfs(boot_3d[[i_rep]], expected, format = "dadi", overwrite = TRUE)
    expect_identical(
      readLines(file.path(
        dadi_dir, paste0("bootstrap_", i_rep), "popA-popB-popC.sfs"
      )),
      readLines(expected)
    )
  }
  # the replicates do not depend on the size of the chunks in which they are
  # projected
  out_dir_chunks <- withr::local_tempdir()
  withr::with_options(list(tidypopgen.sfs_max_chunk_cells = 1), {
    withr::with_seed(20, {
      write_sfs(sim$gt,
        format = "fastsimcoal2",
        out_dir = out_dir_chunks, n_proj = n_proj,
        prefix = "test", n_bootstrap = 3, bootstrap_block_size = 20
      )
    })
  })
  for (f in list.files(out_dir, recursive = TRUE)) {
    expect_identical(
      readLines(file.path(out_dir_chunks, f)),
      readLines(file.path(out_dir, f))
    )
  }
  expect_error(
    write_sfs(sim$gt,
      format = "fastsimcoal2",
      out_dir = out_dir, n_proj = n_proj, prefix = "test",
      n_bootstrap = 3, bootstrap_block_size = 20
    ),
    "already exist"
  )
  expect_error(
    write_sfs(sim$gt,
      out_dir = withr::local_tempdir(), n_proj = n_proj,
      n_bootstrap = 3
    ),
    "`block_size` is needed"
  )
})

test_that("replicates with many blocks do not depend on the number of cores", {
  # 60 blocks of 2 SNPs: the block spectra are buffered (32 at a time) before
  # being added to the replicates, and a single spectrum adds them in
  # parallel
  block_id <- .jackknife_block_ids(show_loci(sim$gt),
    block_size = 2,
    size_unit = "snp"
  )
  boot_1 <- withr::with_seed(30, {
    bootstrap_nwise_pop_sfs(sim$gt,
      n_bootstrap = 4, block_size = 2,
      n_proj = n_proj, type = "array", n_cores = 1
    )
  })
  boot_2 <- withr::with_seed(30, {
    bootstrap_nwise_pop_sfs(sim$gt,
      n_bootstrap = 4, block_size = 2,
      n_proj = n_proj, type = "array", n_cores = 2
    )
  })
  expect_identical(boot_1, boot_2)
  for (i_rep in c(1, 4)) {
    gt_rep <- resampled_gt(
      sim, block_id,
      attr(boot_1, "bootstrap_blocks")[i_rep, ]
    )
    expect_equal(
      boot_1[[i_rep]],
      nwise_pop_sfs(gt_rep, n_proj = n_proj, type = "array"),
      tolerance = 1e-10
    )
  }
})
