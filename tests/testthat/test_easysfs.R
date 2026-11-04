# Compare the site frequency spectra with easySFS
# The vcfs and the easySFS outputs are generated with easySFS (see
# data-raw/create_easysfs_test_data.R for the script and the commit of
# easySFS used). The tests are skipped if the data have not been created.

easysfs_dir <- test_path("testdata", "easysfs")

skip_if_no_easysfs_data <- function() {
  skip_if_not(
    file.exists(file.path(easysfs_dir, "easysfs_version.txt")),
    "easySFS test data not created (see data-raw/create_easysfs_test_data.R)"
  )
}

# read a vcf and its populations file into a grouped gen_tibble, with the
# populations as a factor with the given levels
read_easysfs_gt <- function(vcf, pops_file, levels) {
  gt <- gen_tibble(
    file.path(easysfs_dir, vcf),
    backingfile = tempfile(),
    quiet = TRUE
  )
  pops <- utils::read.table(file.path(easysfs_dir, pops_file),
    col.names = c("id", "population")
  )
  gt$population <- factor(pops$population[match(gt$id, pops$id)],
    levels = levels
  )
  gt %>% dplyr::group_by(population)
}

# compare a tidy spectrum (with the last population changing fastest) to a
# dadi file
expect_tidy_equals_dadi <- function(sfs_tidy, file) {
  dadi_sfs <- read_dadi_sfs(file)
  expect_equal(sfs_tidy$n_sites, dadi_sfs$data, tolerance = 1e-10)
  expect_identical(sfs_tidy$masked, dadi_sfs$mask)
  expect_identical(attr(sfs_tidy, "folded"), dadi_sfs$folded)
  dadi_sfs
}

pops <- c("popA", "popB", "popC")
proj <- c(popA = 10, popB = 8, popC = 6)

test_that("pop_sfs, pairwise_pop_sfs and nwise_pop_sfs match easySFS", {
  skip_if_no_easysfs_data()
  easysfs_gt <- read_easysfs_gt(
    "easysfs_test.vcf", "easysfs_populations.txt",
    levels = pops
  )
  for (folding in c("folded", "unfolded")) {
    folded <- folding == "folded"
    dadi_dir <- file.path(easysfs_dir, folding, "dadi")
    # 1D
    sfs_1d <- pop_sfs(easysfs_gt, n_proj = proj, folded = folded)
    for (pop in pops) {
      dadi_sfs <- expect_tidy_equals_dadi(
        sfs_1d[sfs_1d$group == pop, ],
        file.path(dadi_dir, paste0(pop, "-", proj[[pop]], ".sfs"))
      )
      expect_identical(dadi_sfs$dims, as.integer(proj[[pop]] + 1))
    }
    # 2D
    sfs_2d <- pairwise_pop_sfs(easysfs_gt, n_proj = proj, folded = folded)
    pairs <- utils::combn(pops, 2)
    for (i in seq_len(ncol(pairs))) {
      is_pair <- sfs_2d$population_1 == pairs[1, i] &
        sfs_2d$population_2 == pairs[2, i]
      expect_tidy_equals_dadi(
        sfs_2d[is_pair, ],
        file.path(dadi_dir, paste0(pairs[1, i], "-", pairs[2, i], ".sfs"))
      )
    }
    # multi
    sfs_nd <- nwise_pop_sfs(easysfs_gt, n_proj = proj, folded = folded)
    dadi_sfs <- expect_tidy_equals_dadi(
      sfs_nd,
      file.path(dadi_dir, "popA-popB-popC.sfs")
    )
    expect_identical(dadi_sfs$pop_ids, pops)
  }
})

test_that("the order of the populations matches easySFS --order", {
  skip_if_no_easysfs_data()
  # easySFS was run with --order popC,popA,popB
  ordered_pops <- c("popC", "popA", "popB")
  easysfs_gt <- read_easysfs_gt(
    "easysfs_test.vcf", "easysfs_populations.txt",
    levels = ordered_pops
  )
  dadi_dir <- file.path(easysfs_dir, "ordered", "dadi")
  sfs_2d <- pairwise_pop_sfs(easysfs_gt, n_proj = proj)
  pairs <- utils::combn(ordered_pops, 2)
  for (i in seq_len(ncol(pairs))) {
    is_pair <- sfs_2d$population_1 == pairs[1, i] &
      sfs_2d$population_2 == pairs[2, i]
    expect_tidy_equals_dadi(
      sfs_2d[is_pair, ],
      file.path(dadi_dir, paste0(pairs[1, i], "-", pairs[2, i], ".sfs"))
    )
  }
  sfs_nd <- nwise_pop_sfs(easysfs_gt, n_proj = proj)
  dadi_sfs <- expect_tidy_equals_dadi(
    sfs_nd,
    file.path(dadi_dir, "popC-popA-popB.sfs")
  )
  expect_identical(dadi_sfs$pop_ids, ordered_pops)
})

test_that("pop_sfs_preview matches easySFS --preview", {
  skip_if_no_easysfs_data()
  easysfs_gt <- read_easysfs_gt(
    "easysfs_test.vcf", "easysfs_populations.txt",
    levels = pops
  )
  preview <- pop_sfs_preview(easysfs_gt)
  # easySFS prints the name of each population, followed by a line of
  # (projection, rounded number of segregating sites)
  easysfs_preview <- readLines(file.path(easysfs_dir, "easysfs_preview.txt"))
  for (pop in pops) {
    values_line <- easysfs_preview[which(trimws(easysfs_preview) == pop) + 1]
    values <- regmatches(
      values_line,
      gregexpr("\\(([0-9]+), ([0-9.]+)\\)", values_line)
    )[[1]]
    values <- do.call(rbind, strsplit(gsub("[()]", "", values), ", "))
    this_preview <- preview[preview$group == pop, ]
    expect_identical(this_preview$n_proj, as.integer(values[, 1]))
    # easySFS only prints rounded values: for values half way between two
    # integers, the last bits decide on the rounding, so either is accepted
    n_seg <- this_preview$n_seg_sites
    easysfs_seg <- as.numeric(values[, 2])
    tie <- abs(n_seg - floor(n_seg) - 0.5) < 1e-8
    expect_equal(round(n_seg[!tie]), easysfs_seg[!tie])
    expect_true(all(abs(n_seg[tie] - easysfs_seg[tie]) <= 0.5 + 1e-8))
  }
})

test_that("pop_sfs matches easySFS on haploid data", {
  skip_if_no_easysfs_data()
  haploid_gt <- read_easysfs_gt(
    "easysfs_haploid.vcf", "easysfs_haploid_populations.txt",
    levels = c("popA", "popB")
  )
  hap_proj <- c(popA = 5, popB = 4)
  dadi_dir <- file.path(easysfs_dir, "haploid", "dadi")
  sfs_1d <- pop_sfs(haploid_gt, n_proj = hap_proj)
  for (pop in c("popA", "popB")) {
    expect_tidy_equals_dadi(
      sfs_1d[sfs_1d$group == pop, ],
      file.path(dadi_dir, paste0(pop, "-", hap_proj[[pop]], ".sfs"))
    )
  }
  expect_tidy_equals_dadi(
    pairwise_pop_sfs(haploid_gt, n_proj = hap_proj),
    file.path(dadi_dir, "popA-popB.sfs")
  )
})

# compare all the files written by gt_as_easysfs() and by easySFS (the text,
# and the numbers to a tolerance)
expect_same_files_as_easysfs <- function(our_dir, easysfs_out_dir) {
  easysfs_files <- sort(list.files(easysfs_out_dir, recursive = TRUE))
  expect_identical(sort(list.files(our_dir, recursive = TRUE)), easysfs_files)
  for (f in easysfs_files) {
    ours <- sfs_file_parts(file.path(our_dir, f))
    theirs <- sfs_file_parts(file.path(easysfs_out_dir, f))
    expect_identical(ours$text, theirs$text, label = f)
    expect_equal(ours$numbers, theirs$numbers, tolerance = 1e-10, label = f)
  }
}

test_that("gt_as_easysfs writes the same files as easySFS", {
  skip_if_no_easysfs_data()
  easysfs_gt <- read_easysfs_gt(
    "easysfs_test.vcf", "easysfs_populations.txt",
    levels = pops
  )
  runs <- list(
    folded = list(levels = pops, args = list()),
    unfolded = list(levels = pops, args = list(folded = FALSE)),
    ordered = list(levels = c("popC", "popA", "popB"), args = list()),
    total_length = list(levels = pops, args = list(total_length = 100000))
  )
  for (run in names(runs)) {
    this_gt <- easysfs_gt %>%
      dplyr::ungroup() %>%
      dplyr::mutate(population = factor(population,
        levels = runs[[run]]$levels
      )) %>%
      dplyr::group_by(population)
    out_dir <- withr::local_tempdir()
    do.call(gt_as_easysfs, c(
      list(this_gt, out_dir = out_dir, n_proj = proj, prefix = "test"),
      runs[[run]]$args
    ))
    expect_same_files_as_easysfs(out_dir, file.path(easysfs_dir, run))
  }
  # haploid data
  haploid_gt <- read_easysfs_gt(
    "easysfs_haploid.vcf", "easysfs_haploid_populations.txt",
    levels = c("popA", "popB")
  )
  out_dir <- withr::local_tempdir()
  gt_as_easysfs(haploid_gt,
    out_dir = out_dir, n_proj = c(popA = 5, popB = 4),
    prefix = "test"
  )
  expect_same_files_as_easysfs(out_dir, file.path(easysfs_dir, "haploid"))
})
