raw_path_pop_a <- system.file("extdata/pop_a.raw", package = "tidypopgen")
map_path_pop_a <- system.file("extdata/pop_a.map", package = "tidypopgen")
pop_a_gt <- read_plink_raw(file = raw_path_pop_a, map_file = map_path_pop_a, quiet = TRUE)
out_file <- tempfile()

testthat::test_that("we can export to plink raw",{
  # delete any files if they already exists (unlikely)
  unlink(paste0(out_file,"*"))
  expect_true(gt_as_plink(pop_a_gt, file = out_file, chunk_size = 2))
  out_file_raw <- paste0(out_file,".raw")
  expect_true(file.exists(out_file_raw))
  # check that the file that we generated is identical to the original
  expect_true(tools::md5sum(out_file_raw)==tools::md5sum(raw_path_pop_a))

  out_file_map <- paste0(out_file,".map")
  expect_true(file.exists(out_file_map))
  new_pop_a_gt <- read_plink_raw(file = out_file_raw,
                                 map_file = out_file_map, quiet = TRUE)
  expect_identical(pop_a_gt, new_pop_a_gt)
})

testthat::test_that("we can export to plink ped",{
  # delete any files if they already exists
  unlink(paste0(out_file,"*"))
  expect_true(gt_as_plink(pop_a_gt, file = out_file, plink_format = "ped", chunk_size = 2))
  out_file_raw <- paste0(out_file,".ped")
  expect_true(file.exists(out_file_raw))
})
