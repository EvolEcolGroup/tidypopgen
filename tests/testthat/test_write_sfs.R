# Tests for write_sfs(), for sfs objects and gen_tibbles

# a small unfolded 2D spectrum (populations a and b, with 1 and 2 copies)
sfs_2d <- new_sfs(
  array(c(0, 1.5, 2, 3, 4.25, 5), dim = c(2, 3)),
  mask = array(c(TRUE, rep(FALSE, 4), TRUE), dim = c(2, 3)),
  pop_ids = c("a", "b"),
  folded = FALSE,
  n_loci = 15.75,
  n_loci_dropped = 0
)
# a small folded 1D spectrum
sfs_1d <- new_sfs(
  array(c(1, 1 / 3, 2), dim = 3),
  mask = array(c(TRUE, FALSE, TRUE), dim = 3),
  pop_ids = "a",
  folded = TRUE,
  n_loci = 10 / 3,
  n_loci_dropped = 1
)

test_that("write_sfs writes dadi files", {
  file <- withr::local_tempfile(fileext = ".sfs")
  expect_identical(write_sfs(sfs_2d, file), file)
  # the last population changes fastest
  expect_identical(readLines(file), c(
    "2 3 unfolded \"a\" \"b\"",
    "0 2 4.25 1.5 3 5",
    "1 0 0 0 0 1"
  ))
  # the first population changes fastest, with the populations reversed
  write_sfs(sfs_2d, file, fastest = "first", overwrite = TRUE)
  expect_identical(readLines(file), c(
    "3 2 unfolded \"b\" \"a\"",
    "0 1.5 2 3 4.25 5",
    "1 0 0 0 0 1"
  ))
  # numbers are written as numpy.savetxt(fmt = "%.16g")
  write_sfs(sfs_1d, file, overwrite = TRUE)
  expect_identical(readLines(file), c(
    "3 folded \"a\"",
    "1 0.3333333333333333 2",
    "1 0 1"
  ))
})

test_that("write_sfs writes fastsimcoal2 files", {
  file <- withr::local_tempfile(fileext = ".obs")
  write_sfs(sfs_1d, file, format = "fastsimcoal2")
  expect_identical(readLines(file), c(
    "1 observation",
    "d0_0\td0_1\td0_2",
    "1 0.3333333333333333 2",
    ""
  ))
  # 2D: the second population in the rows, the first in the columns
  write_sfs(sfs_2d, file, format = "fastsimcoal2", overwrite = TRUE)
  expect_identical(readLines(file), c(
    "1 observation",
    "\td0_0\td0_1",
    "d1_0\t0 1.5",
    "d1_1\t2 3",
    "d1_2\t4.25 5"
  ))
  # multidimensional
  write_sfs(sfs_2d, file,
    format = "fastsimcoal2_multidimensional",
    overwrite = TRUE
  )
  expect_identical(readLines(file), c(
    "1 observations. No. of demes and sample sizes are on next line.",
    "2\t1 2",
    "0 2 4.25 1.5 3 5",
    ""
  ))
  write_sfs(sfs_2d, file,
    format = "fastsimcoal2_multidimensional", fastest = "first",
    overwrite = TRUE
  )
  expect_identical(readLines(file)[2:3], c("2\t1 2", "0 1.5 2 3 4.25 5"))
  # a spectrum of three populations is written in the multidimensional format
  sfs_3d <- load_example_gt("grouped_gen_tbl") %>%
    nwise_pop_sfs(n_proj = 2, type = "array")
  write_sfs(sfs_3d, file, format = "fastsimcoal2", overwrite = TRUE)
  fsc_lines <- readLines(file)
  expect_identical(fsc_lines[2], "3\t2 2 2")
  expect_equal(
    as.numeric(strsplit(fsc_lines[3], " ")[[1]]),
    tidy(sfs_3d)$n_sites
  )
})

test_that("write_sfs sets the monomorphic sites and rounds", {
  file <- withr::local_tempfile(fileext = ".sfs")
  write_sfs(sfs_2d, file, total_length = 100)
  # the first cell is the total length minus the sites in the other cells
  expect_identical(readLines(file)[2], "84.25 2 4.25 1.5 3 5")
  write_sfs(sfs_2d, file, total_length = 100, round = TRUE, overwrite = TRUE)
  expect_identical(readLines(file)[2], "84 2 4 2 3 5")
  expect_error(
    write_sfs(sfs_2d, file, total_length = 15, overwrite = TRUE),
    "should be larger than the number of sites"
  )
  expect_error(
    write_sfs(sfs_2d, file, total_length = "a", overwrite = TRUE),
    "single number"
  )
})

test_that("write_sfs checks its arguments", {
  file <- withr::local_tempfile(fileext = ".sfs")
  write_sfs(sfs_2d, file)
  expect_error(write_sfs(sfs_2d, file), "already exists")
  # formats are named in full
  expect_error(
    write_sfs(sfs_2d, file, format = "fsc", overwrite = TRUE),
    "should be one of"
  )
  expect_error(write_sfs(unclass(sfs_2d), file), "should be an sfs object")
  expect_error(
    write_sfs(sfs_2d, file, round = NA, overwrite = TRUE),
    "`round` should be TRUE or FALSE"
  )
})

test_that("write_sfs files can be read back", {
  sim <- ref_sfs_simulate(c(popA = 5, popB = 4, popC = 3),
    n_loci = 60,
    seed = 31
  )
  sfs <- nwise_pop_sfs(sim$gt, n_proj = c(6, 4, 3), type = "array")
  file <- withr::local_tempfile(fileext = ".sfs")
  write_sfs(sfs, file)
  dadi_sfs <- read_dadi_sfs(file)
  expect_identical(dadi_sfs$dims, c(7L, 5L, 4L))
  expect_identical(dadi_sfs$pop_ids, c("popA", "popB", "popC"))
  expect_true(dadi_sfs$folded)
  sfs_tidy <- tidy(sfs)
  # %.16g keeps 16 significant digits
  expect_equal(dadi_sfs$data, sfs_tidy$n_sites, tolerance = 1e-15)
  expect_identical(dadi_sfs$mask, sfs_tidy$masked)
})

test_that("write_sfs writes the files of easySFS for a gen_tibble", {
  sim <- ref_sfs_simulate(c(popA = 5, popB = 4, popC = 3),
    n_loci = 60,
    seed = 32
  )
  dadi_dir <- withr::local_tempdir()
  fsc_dir <- withr::local_tempdir()
  n_proj <- c(popA = 6, popB = 4, popC = 3)
  dadi_files <- write_sfs(sim$gt, out_dir = dadi_dir, n_proj = n_proj)
  fsc_files <- write_sfs(sim$gt,
    format = "fastsimcoal2",
    out_dir = fsc_dir, n_proj = n_proj,
    prefix = "test"
  )
  expect_identical(sort(basename(dadi_files)), sort(c(
    "popA-6.sfs", "popB-4.sfs", "popC-3.sfs", "popA-popB.sfs",
    "popA-popC.sfs", "popB-popC.sfs", "popA-popB-popC.sfs"
  )))
  expect_identical(sort(basename(fsc_files)), sort(c(
    "popA_MAFpop0.obs", "popB_MAFpop0.obs", "popC_MAFpop0.obs",
    "test_jointMAFpop1_0.obs", "test_jointMAFpop2_0.obs",
    "test_jointMAFpop2_1.obs", "test_MSFS.obs"
  )))
  # the files are written directly in out_dir
  expect_identical(dirname(dadi_files), rep(dadi_dir, 7))
  expect_identical(sort(list.files(dadi_dir)), sort(basename(dadi_files)))
  expect_identical(sort(list.files(fsc_dir)), sort(basename(fsc_files)))
  # the files are those written by write_sfs() for each spectrum
  expect_same_file <- function(sfs, format, file, ...) {
    tmp <- withr::local_tempfile()
    write_sfs(sfs, tmp, format = format, ...)
    expect_identical(readLines(file), readLines(tmp))
  }
  sfs_1d <- pop_sfs(sim$gt, n_proj = n_proj, type = "array")
  sfs_2d <- pairwise_pop_sfs(sim$gt, n_proj = n_proj, type = "array")
  sfs_3d <- nwise_pop_sfs(sim$gt, n_proj = n_proj, type = "array")
  expect_same_file(sfs_1d$popB, "dadi", file.path(dadi_dir, "popB-4.sfs"))
  expect_same_file(
    sfs_1d$popB, "fastsimcoal2",
    file.path(fsc_dir, "popB_MAFpop0.obs")
  )
  expect_same_file(
    sfs_2d$popA.popC, "dadi",
    file.path(dadi_dir, "popA-popC.sfs")
  )
  expect_same_file(
    sfs_3d, "dadi",
    file.path(dadi_dir, "popA-popB-popC.sfs")
  )
  expect_same_file(
    sfs_3d, "fastsimcoal2_multidimensional",
    file.path(fsc_dir, "test_MSFS.obs")
  )
  # the joint spectrum of popA (0) and popC (2) has popC in the rows
  fsc_2d <- readLines(file.path(fsc_dir, "test_jointMAFpop2_0.obs"))
  col_labels <- paste0("d0_", 0:6, collapse = "\t")
  expect_identical(fsc_2d[2], paste0("\t", col_labels))
  expect_identical(
    substr(fsc_2d[3:6], 1, 5),
    c("d2_0\t", "d2_1\t", "d2_2\t", "d2_3\t")
  )
  expect_equal(
    as.numeric(strsplit(sub("^d2_1\t", "", fsc_2d[4]), " ")[[1]]),
    as.vector(sfs_2d$popA.popC[, 2])
  )

  # files are not overwritten
  expect_error(
    write_sfs(sim$gt, out_dir = dadi_dir, n_proj = n_proj),
    "already exist"
  )
  expect_error(
    write_sfs(sim$gt,
      format = "fastsimcoal2", out_dir = fsc_dir, n_proj = n_proj,
      prefix = "test"
    ),
    "already exist"
  )
  # unfolded spectra, the first population changing fastest, the number of
  # monomorphic sites, and rounding
  dadi_files <- write_sfs(sim$gt,
    out_dir = dadi_dir, n_proj = n_proj, folded = FALSE,
    fastest = "first", total_length = 1000, round = TRUE, overwrite = TRUE
  )
  fsc_files <- write_sfs(sim$gt,
    format = "fastsimcoal2",
    out_dir = fsc_dir, n_proj = n_proj, prefix = "test", folded = FALSE,
    fastest = "first", total_length = 1000, round = TRUE, overwrite = TRUE
  )
  unfolded_files <- c("popA_DAFpop0.obs", "test_jointDAFpop1_0.obs")
  expect_true(all(unfolded_files %in% basename(fsc_files)))
  sfs_unfolded <- nwise_pop_sfs(sim$gt,
    n_proj = n_proj, folded = FALSE,
    type = "array"
  )
  expect_same_file(sfs_unfolded, "dadi",
    file.path(dadi_dir, "popA-popB-popC.sfs"),
    fastest = "first", total_length = 1000, round = TRUE
  )
  expect_same_file(sfs_unfolded, "fastsimcoal2_multidimensional",
    file.path(fsc_dir, "test_MSFS.obs"),
    fastest = "first", total_length = 1000, round = TRUE
  )
})

test_that("write_sfs uses the backing files of a gen_tibble by default", {
  test_gt <- load_example_gt("grouped_gen_tbl")
  backing_name <- bigstatsr::sub_bk(attr(test_gt$genotypes, "fbm")$backingfile)
  dadi_dir <- paste0(backing_name, "_dadi")
  fsc_dir <- paste0(backing_name, "_fastsimcoal2")
  on.exit(unlink(c(dadi_dir, fsc_dir), recursive = TRUE), add = TRUE)
  expect_message(
    dadi_files <- write_sfs(test_gt),
    "maximises its number of segregating sites"
  )
  expect_true(all(startsWith(dadi_files, dadi_dir)))
  fsc_files <- suppressMessages(write_sfs(test_gt, format = "fastsimcoal2"))
  expect_true(all(startsWith(fsc_files, fsc_dir)))
  expect_true(file.exists(file.path(
    fsc_dir,
    paste0(basename(backing_name), "_MSFS.obs")
  )))
})

test_that("write_sfs works with a gen_tibble with a single group", {
  test_gt <- load_example_gt("grouped_gen_tbl") %>%
    dplyr::filter(population == "pop1")
  dadi_files <- write_sfs(test_gt,
    out_dir = withr::local_tempdir(),
    n_proj = 4
  )
  expect_identical(sort(basename(dadi_files)), c("pop1-4.sfs", "pop1.sfs"))
  fsc_files <- write_sfs(test_gt,
    format = "fastsimcoal2",
    out_dir = withr::local_tempdir(),
    n_proj = 4, prefix = "test"
  )
  expect_identical(
    sort(basename(fsc_files)),
    c("pop1_MAFpop0.obs", "test_MSFS.obs")
  )
})

test_that("write_sfs checks its arguments for a gen_tibble", {
  test_gt <- load_example_gt("grouped_gen_tbl")
  out_dir <- withr::local_tempdir()
  expect_error(
    write_sfs(dplyr::ungroup(test_gt), out_dir = out_dir),
    "grouped gen_tibble"
  )
  expect_error(
    write_sfs(test_gt,
      format = "fastsimcoal2", out_dir = out_dir, n_proj = 2,
      max_cells = 10
    ),
    "27 cells"
  )
  test_gt_slash <- test_gt %>%
    dplyr::ungroup() %>%
    dplyr::mutate(population = paste0(population, "/x")) %>%
    dplyr::group_by(population)
  expect_error(
    write_sfs(test_gt_slash, out_dir = out_dir, n_proj = 2),
    "cannot contain"
  )
  expect_error(
    write_sfs(test_gt,
      format = "fastsimcoal2", out_dir = out_dir, n_proj = 2,
      fastest = "middle"
    ),
    "should be one of"
  )
})
