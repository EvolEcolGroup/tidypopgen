# limit number of threads for tests
data.table::setDTthreads(2)
if (rlang::is_installed("RhpcBLASctl")) {
  RhpcBLASctl::blas_set_num_threads(2)
  RhpcBLASctl::omp_set_num_threads(2)
}

test_that("cbind correctly binds", {
  test_genotypes <- rbind(
    c(1, 1, 0, 1, 1, 0),
    c(2, 1, 0, NA, 0, 0),
    c(2, NA, 0, 0, 1, 1),
    c(1, 0, 0, 1, 0, 0),
    c(1, 2, 0, 1, 2, 1),
    c(0, 0, 0, 0, NA, 1),
    c(0, 1, 1, 0, 1, NA)
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
  # new data.frame to cbind
  new_df <- data.frame(
    id2 = c("a", "b", "c", "d", "e", "f", "g"),
    population = c(
      "region1", "region1", "region2", "region2",
      "region1", "region3", "region3"
    ),
    age = c(1, 2, 3, 4, 5, 6, 7)
  )
  # population is in both, so both columns are renamed (as in bind_cols())
  expect_message(
    test_combined_gt <- cbind(test_gt, new_df),
    "New names"
  )
  expect_true(inherits(test_combined_gt, "gen_tbl"))
  expect_identical(
    names(test_combined_gt),
    c(
      "id", "population...2", "genotypes", "id2", "population...5", "age"
    )
  )
  # check that the loci info is there
  expect_equal(nrow(show_loci(test_combined_gt)), nrow(test_loci))

  # a grouped gen_tibble keeps its grouping, unless its grouping variable is
  # renamed; the names are those of bind_cols()
  test_grouped <- test_gt %>% group_by(population)
  expect_message(
    test_combined_gt <- cbind(test_grouped, new_df),
    "New names"
  )
  expect_false(dplyr::is_grouped_df(test_combined_gt))
  expect_true(inherits(test_combined_gt, "gen_tbl"))
  expect_identical(
    names(test_combined_gt),
    names(suppressMessages(dplyr::bind_cols(
      dplyr::group_by(tibble::as_tibble(test_indiv_meta), population),
      genotypes = 1,
      new_df
    )))
  )
  test_combined_gt <- cbind(test_grouped, new_df["age"])
  expect_identical(dplyr::group_vars(test_combined_gt), "population")
  expect_true(inherits(test_combined_gt, "grouped_gen_tbl"))
  test_combined_gt <- suppressMessages(
    cbind(test_grouped, data.frame(age = 1:7, age = 7:1, check.names = FALSE))
  )
  expect_identical(dplyr::group_vars(test_combined_gt), "population")
  expect_identical(
    names(test_combined_gt),
    c("id", "population", "genotypes", "age...4", "age...5")
  )

  # the id and genotypes can not be duplicated
  expect_error(
    cbind(test_gt, data.frame(id = test_gt$id)),
    "already has a column 'id'"
  )
  expect_error(
    cbind(test_grouped, genotypes = 1:7),
    "already has a column 'genotypes'"
  )
})
