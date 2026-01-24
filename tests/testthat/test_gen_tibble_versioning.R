test_that("versioning updates correctly through gt_order_loci", {
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c"),
    population = c("pop1", "pop1", "pop2")
  )
  test_genotypes <- rbind(
    c(1, 1, 0, 1, 1, 0),
    c(2, 1, 0, 0, 0, 0),
    c(2, 2, 0, 0, 1, 1)
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
    quiet = TRUE,
    backingfile = tempfile()
  )

  # get the gt filenames
  files <- gt_get_file_names(test_gt)

  # do the files exist?
  expect_true(file.exists(files[1]))
  expect_true(file.exists(files[2]))

  # remove extension
  file <- gsub(".bk", "", files[2], fixed = TRUE)

  # does the directory still exist?
  tempfile_dir <- dirname(file)
  expect_true(dir.exists(tempfile_dir))

  # create gt using the same backingfile name
  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE,
    backingfile = file
  )

  # get new file names
  new_files <- gt_get_file_names(test_gt)

  # mess with the loci table
  test_gt <- test_gt %>% select_loci(c(2, 4, 5, 1, 6))
  test_gt <- gt_order_loci(test_gt, use_current_table = FALSE, quiet = TRUE)
  gt_save(test_gt, quiet = TRUE)

  # check file names
  expect_equal(gt_get_file_names(test_gt)[1], paste0(file, "_v3.rds"))
  expect_equal(gt_get_file_names(test_gt)[2], paste0(file, "_v3.bk"))
})


test_that("versioning if .bk already exists", {
  test_indiv_meta <- data.frame(
    id = c("a", "b", "c"),
    population = c("pop1", "pop1", "pop2")
  )
  test_genotypes <- rbind(
    c(1, 1, 0, 1, 1, 0),
    c(2, 1, 0, 0, 0, 0),
    c(2, 2, 0, 0, 1, 1)
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
    quiet = TRUE,
    backingfile = tempfile()
  )

  # get the gt filenames
  files <- gt_get_file_names(test_gt)

  # do the files exist?
  expect_true(file.exists(files[1]))
  expect_true(file.exists(files[2]))

  # remove the .rds
  expect_true(file.remove(files[1]))

  # remove extension
  file <- gsub(".bk", "", files[2], fixed = TRUE)

  # does the directory still exist?
  tempfile_dir <- dirname(file)
  expect_true(dir.exists(tempfile_dir))

  # create gt using the same backingfile name
  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE,
    backingfile = file
  )

  # get new file names
  new_files <- gt_get_file_names(test_gt)

  # new_files has the same name as original file, plus a version extension
  expect_equal(new_files[2], paste0(file, "_v2.bk"))

  # repeating the process creates another version
  expect_true(file.remove(new_files[1]))
  expect_false(file.exists(new_files[1]))
  expect_true(file.exists(new_files[2]))

  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE,
    backingfile = file
  )

  new_version_files <- gt_get_file_names(test_gt)

  expect_equal(new_version_files[2], paste0(file, "_v3.bk"))
})
