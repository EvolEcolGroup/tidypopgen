# create file
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
  quiet = TRUE
)

# this also tests show_genotypes and show_loci
test_that("save and load gt", {
  expect_true(inherits(test_gt, "gen_tbl"))
  # now save the tibble
  all_file_names <- gt_save(test_gt, quiet = TRUE)
  # check that the new file exists
  expect_true(file.exists(all_file_names[1]))
  new_test_gt <- gt_load(all_file_names[1])

  # check that we preserved the genotypes
  expect_true(all(show_genotypes(new_test_gt$genotypes) == test_genotypes))
  # check that we preserved the loci
  expect_identical(
    show_loci(new_test_gt$genotypes) %>%
      select(-big_index, -chromosome),
    as_tibble(test_loci[, c(
      "name",
      "position",
      "genetic_dist",
      "allele_ref",
      "allele_alt"
    )])
  ) # nolint
  expect_identical(
    show_loci(test_gt$genotypes) %>%
      select(c(-big_index, -chromosome)),
    as_tibble(test_loci[, c(
      "name",
      "position",
      "genetic_dist",
      "allele_ref",
      "allele_alt"
    )])
  ) # nolint

  # now remove the tibble
  rm(new_test_gt)
  # now move the backing files
  new_dir <- file.path(dirname(all_file_names[1]), "test")
  dir.create(new_dir)
  expect_true(file.copy(
    from = all_file_names[2],
    to = file.path(new_dir, basename(all_file_names[2]))
  ))
  expect_true(file.copy(
    from = all_file_names[3],
    to = file.path(new_dir, basename(all_file_names[3]))
  ))
  expect_true(file.remove(all_file_names[2]))

  # expect_true(file.remove(all_file_names[3])) #nolint
  # TODO the above test fails on Windows,
  # needs a fix after response to bigstatsr issue

  # loading should fail
  expect_error(
    new_test_gt2 <- gt_load(all_file_names[1]),
    ".rds does not exist"
  )
  # this should now work:
  new_test_gt2 <-
    gt_load(
      all_file_names[1],
      reattach_to = file.path(new_dir, basename(all_file_names[2]))
    )
  # verify that we have all the info
  # check that we preserved the genotypes
  expect_true(all(show_genotypes(new_test_gt2$genotypes) == test_genotypes))
  # check that we preserved the loci
  expect_identical(
    show_loci(new_test_gt2$genotypes) %>%
      select(-big_index, -chromosome),
    as_tibble(test_loci[, c(
      "name",
      "position",
      "genetic_dist",
      "allele_ref",
      "allele_alt"
    )])
  ) # nolint
  expect_identical(
    show_loci(test_gt$genotypes) %>% select(-big_index, -chromosome),
    as_tibble(test_loci[, c(
      "name",
      "position",
      "genetic_dist",
      "allele_ref",
      "allele_alt"
    )])
  ) # nolint
})

test_that("error if saving a non gen_tibble object", {
  expect_error(
    gt_save(
      indiv_meta,
      "x should be a gen_tibble"
    )
  )
})

test_that("old bigsnp gen_tibbles can be loaded", {
  # TODO find a way to run this test on CI
  # copy the old gen_tibble to a temp location  #nolint start
  temp_dir <- tempdir()
  file.copy(system.file("extdata/bigsnp_obj_gt/old_gen_tibble.gt",
    package = "tidypopgen"
  ), to = file.path(temp_dir, "old_gen_tibble.gt"))
  file.copy(system.file("extdata/bigsnp_obj_gt/old_gen_tibble.bk",
    package = "tidypopgen"
  ), to = file.path(temp_dir, "old_gen_tibble.bk"))
  file.copy(system.file("extdata/bigsnp_obj_gt/old_gen_tibble.rds",
    package = "tidypopgen"
  ), to = file.path(temp_dir, "old_gen_tibble.rds"))

  gt_path <- paste0(temp_dir, "/old_gen_tibble.gt")
  rds_path <- paste0(temp_dir, "/old_gen_tibble.rds")
  bk_path <- paste0(temp_dir, "/old_gen_tibble.bk")

  # The test works without editing attributes, but it updates the .rds in
  # extdata each time

  # To fix, we want to update attributes

  # load gt using readRDS
  gt <- readRDS(gt_path)

  # change paths in gt attributes
  attributes(gt$genotypes)$bigsnp_file <- rds_path
  attributes(gt$genotypes)$bigsnp <- bigsnpr::snp_attach(rds_path)


  # save
  new_rds <- saveRDS(gt, gt_path)

  # reload with gt_load
  expect_message(
    test_gt <-
      gt_load(gt_path), "your gen_tibble was in an old format"
  )
  example_gt <- load_example_gt("gen_tbl")
  expect_true(all.equal(test_gt, example_gt, check.attributes = FALSE))
  # check we have fbm_ploidy
  expect_equal(indiv_ploidy(example_gt), indiv_ploidy(test_gt))
  expect_true(length(indiv_ploidy(test_gt)) == nrow(test_gt))
  expect_true(length(attr(test_gt$genotypes, "fbm_ploidy")) == nrow(test_gt))
  # check it works for basic calculations
  global_stats <- test_gt %>% pop_global_stats()
  global_stats_example <- example_gt %>% pop_global_stats()
  expect_true(all.equal(global_stats,
    global_stats_example,
    check.attributes = FALSE
  )) # nolint end
})
