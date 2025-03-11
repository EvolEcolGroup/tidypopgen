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
    show_loci(new_test_gt$genotypes) %>% select(-big_index, -chr_int),
    as_tibble(test_loci)
  ) # nolint
  expect_identical(
    show_loci(test_gt$genotypes) %>% select(c(-big_index, -chr_int)),
    as_tibble(test_loci)
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
    show_loci(new_test_gt2$genotypes) %>% select(-big_index, -chr_int),
    as_tibble(test_loci)
  ) # nolint
  expect_identical(
    show_loci(test_gt$genotypes) %>% select(-big_index, -chr_int),
    as_tibble(test_loci)
  ) # nolint
})
