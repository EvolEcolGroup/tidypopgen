test_that("loci_hwe produces the same output as plink --hardy midp", {
  # Create gentibble for our data
  bed_path <- system.file(
    "extdata/related/families.bed",
    package = "tidypopgen"
  )
  families_bigsnp_path <- bigsnpr::snp_readBed(
    bed_path,
    backingfile = tempfile()
  )
  families <- gen_tibble(
    families_bigsnp_path,
    quiet = TRUE,
    valid_alleles = c("1", "2")
  )

  # Read in plink results (first 10 snps)
  plink_hwe <- read.table(
    system.file(
      "extdata/related/families_hwe_midp.hwe",
      package = "tidypopgen"
    ),
    header = TRUE
  )

  # Calculate hwe in tidypopgen
  tidy_hwe <- families %>%
    select_loci(c(1:10)) %>%
    loci_hwe(mid_p = TRUE)

  # Compare
  result <- all.equal(plink_hwe$P, tidy_hwe, tolerance = 0.0001)

  # Check results are the same to 4 decimals
  expect_true(result)
})


test_that("loci_hwe mid_p = FALSE produces the same output as plink --hardy ", {
  # Create gentibble for our data
  bed_path <- system.file(
    "extdata/related/families.bed",
    package = "tidypopgen"
  )
  families_bigsnp_path <- bigsnpr::snp_readBed(
    bed_path,
    backingfile = tempfile()
  )
  families <- gen_tibble(
    families_bigsnp_path,
    quiet = TRUE,
    valid_alleles = c("1", "2")
  )

  # Read in plink results (first 10 snps)
  plink_hwe <- read.table(
    system.file("extdata/related/families_hwe.hwe", package = "tidypopgen"),
    header = TRUE
  )

  # Calculate hwe in tidypopgen
  tidy_hwe <- families %>%
    select_loci(c(1:10)) %>%
    loci_hwe(mid_p = FALSE)

  # Compare
  result <- all.equal(plink_hwe$P, tidy_hwe, tolerance = 0.0001)

  # Check results are the same to 4 decimals
  expect_true(result)
})


test_that("loci_hwe for grouped tibble gives the correct result", {
  test_gt <- load_example_gt("grouped_gen_tbl")

  # compute using .grouped_df method
  list <- loci_hwe(test_gt, type = "list")
  matrix <- loci_hwe(test_gt, type = "matrix")
  tidy <- loci_hwe(test_gt, type = "tidy")

  expect_equal(list[1][[1]], as.vector(matrix[, "pop1"]))
  expect_equal(rownames(matrix), show_loci(test_gt)$name)
  expect_equal(colnames(matrix), group_keys(test_gt)$population)
  tidy_pop1 <- tidy %>%
    filter(group == "pop1") %>%
    select(value)
  expect_equal(list[1][[1]], tidy_pop1$value)

  # subset
  test_gt_subset <- test_gt %>% select_loci(c(1, 2, 3, 4))
  list <- loci_hwe(test_gt_subset, type = "list")
  matrix <- loci_hwe(test_gt_subset, type = "matrix")
  tidy <- loci_hwe(test_gt_subset, type = "tidy")
  expect_equal(list[1][[1]], as.vector(matrix[, "pop1"]))
  expect_equal(rownames(matrix), show_loci(test_gt_subset)$name)
  expect_equal(colnames(matrix), group_keys(test_gt_subset)$population)
  tidy_pop1 <- tidy %>%
    filter(group == "pop1") %>%
    select(value)
  expect_equal(list[1][[1]], tidy_pop1$value)

  # compute using group map
  loci_hwe_map <- test_gt %>% group_map(.f = ~ loci_hwe(.x))
  # use fast cpp code (limit cores to 2)
  loci_hwe_grp <- test_gt %>% loci_hwe(n_cores = 2)
  loci_hwe_grp_pop1 <- loci_hwe_grp %>%
    filter(group == "pop1") %>%
    select(value)
  expect_true(all.equal(loci_hwe_map[1][[1]], loci_hwe_grp_pop1$value))

  # and now with reframe
  loci_hwe_reframe <- test_gt %>%
    reframe(loci_hwe = loci_hwe(genotypes))
  loci_hwe_direct <- test_gt %>%
    loci_hwe(n_cores = 2) %>%
    arrange(group)
  expect_equal(loci_hwe_reframe$loci_hwe, loci_hwe_direct$value)
  # check that the direct method can take a column genotypes
  loci_hwe_direct2 <- test_gt %>%
    loci_hwe(genotypes) %>%
    arrange(group)
  expect_equal(loci_hwe_reframe$loci_hwe, loci_hwe_direct2$value)

  # test a second grouping variable
  test_gt$region <- c("a", "a", "b", "b", "a", "b", "b")
  test_gt <- test_gt %>% group_by(population, region)
  expect_error(
    test_gt %>% loci_hwe(),
    "only works with one grouping variable"
  )
})
