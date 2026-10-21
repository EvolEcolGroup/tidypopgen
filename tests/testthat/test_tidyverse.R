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
test_that("inheritance of group tibbles", {
  # test operations on standard gen_tibbles
  expect_true(inherits(test_gt, "gen_tbl"))
  # add column
  test_gt <- test_gt %>% mutate(x = c(1, 2, 3))
  expect_true(inherits(test_gt, "gen_tbl"))
  # remove column
  test_sub_gt <- test_gt %>% select(-x)
  expect_true(inherits(test_sub_gt, "gen_tbl"))
  # remove genotype column
  test_sub_gt <- test_gt %>% select(-genotypes)
  expect_false(inherits(test_sub_gt, "gen_tbl"))
  # filter
  test_sub_gt <- test_gt %>% filter(population == "pop1")
  expect_true(inherits(test_sub_gt, "gen_tbl"))

  # now group it
  test_group_gt <- test_gt %>% group_by(population)
  expect_true(inherits(test_group_gt, "gen_tbl"))
  expect_true(inherits(test_group_gt, "grouped_gen_tbl"))
  # slice rows
  test_group_gt <- test_group_gt %>% filter(population == "pop1")
  expect_true(inherits(test_group_gt, "gen_tbl"))
  expect_true(inherits(test_group_gt, "grouped_gen_tbl"))
  # add columns
  test_group_gt <- test_group_gt %>% mutate(n = c(1, 2))
  expect_true(inherits(test_group_gt, "gen_tbl"))
  expect_true(inherits(test_group_gt, "grouped_gen_tbl"))

  test_group <- test_group_gt %>% select(-genotypes)
  expect_false(inherits(test_group, "gen_tbl"))
  expect_false(inherits(test_group, "grouped_gen_tbl"))

  test_ungroup_gt <- ungroup(test_group_gt)
  expect_true(inherits(test_ungroup_gt, "gen_tbl"))
})


test_that("sf and grouped methods work", {
  # Load tibbles
  grouped_gen_tbl <- load_example_gt("grouped_gen_tbl")
  grouped_gen_tbl_sf <- load_example_gt("grouped_gen_tbl_sf")
  gen_tbl_sf <- load_example_gt("gen_tbl_sf")

  # Filter
  filtered_group <- grouped_gen_tbl %>% filter(id %in% c("a", "c"))
  expect_equal(class(grouped_gen_tbl), class(filtered_group))
  filtered_group_sf <- grouped_gen_tbl_sf %>% filter(id %in% c("a", "c"))
  expect_equal(class(grouped_gen_tbl_sf), class(filtered_group_sf))
  filtered_sf <- gen_tbl_sf %>% filter(id %in% c("a", "c"))
  expect_equal(class(gen_tbl_sf), class(filtered_sf))

  # Arrange
  arranged_group <- grouped_gen_tbl %>% arrange(id)
  expect_equal(class(grouped_gen_tbl), class(arranged_group))
  arranged_group_sf <- grouped_gen_tbl_sf %>% arrange(id)
  expect_equal(class(grouped_gen_tbl_sf), class(arranged_group_sf))
  arranged_sf <- gen_tbl_sf %>% arrange(id)
  expect_equal(class(gen_tbl_sf), class(arranged_sf))

  # Mutate
  mutated_group <- grouped_gen_tbl %>% mutate(region = "East")
  expect_equal(class(grouped_gen_tbl), class(mutated_group))
  mutated_group_sf <- grouped_gen_tbl_sf %>% mutate(region = "East")
  expect_equal(class(grouped_gen_tbl_sf), class(mutated_group_sf))
  mutated_sf <- gen_tbl_sf %>% mutate(region = "East")
  expect_equal(class(gen_tbl_sf), class(mutated_sf))

  # Cbind
  df <- data.frame(region = c("A", "A", "B", "B", "A", "B", "B"))
  gen_tbl_sf_cbind <- cbind(gen_tbl_sf, df)
  expect_equal(class(gen_tbl_sf), class(gen_tbl_sf_cbind))

  # Assignment "$<-"
  class_before <- class(grouped_gen_tbl)
  grouped_gen_tbl$region <- "East"
  expect_equal(class(grouped_gen_tbl), class_before)
  class_before_sf <- class(grouped_gen_tbl_sf)
  grouped_gen_tbl_sf$region <- "East"
  expect_equal(class(grouped_gen_tbl_sf), class_before_sf)
  class_before_sf_un <- class(gen_tbl_sf)
  gen_tbl_sf$region <- "East"
  expect_equal(class(gen_tbl_sf), class_before_sf_un)
})

test_that("$<- on grouping variables preserves grouped gen_tibble classes", {
  grouped_gen_tbl <- load_example_gt("grouped_gen_tbl")
  grouped_gen_tbl_sf <- load_example_gt("grouped_gen_tbl_sf")

  # modify the grouping variable
  class_before <- class(grouped_gen_tbl)
  new_pops <- rev(grouped_gen_tbl$population)
  grouped_gen_tbl$population <- new_pops
  expect_equal(class(grouped_gen_tbl), class_before)
  # groups are recomputed
  expect_equal(group_vars(grouped_gen_tbl), "population")
  expect_equal(
    group_indices(grouped_gen_tbl),
    as.integer(factor(new_pops))
  )
  # and we can still filter
  filtered_group <- grouped_gen_tbl %>% filter(population == "pop1")
  expect_equal(class(filtered_group), class_before)

  # same with sf
  class_before_sf <- class(grouped_gen_tbl_sf)
  grouped_gen_tbl_sf$population <- rev(grouped_gen_tbl_sf$population)
  expect_equal(class(grouped_gen_tbl_sf), class_before_sf)

  # removing the grouping variable returns an ungrouped gen_tibble
  grouped_gen_tbl$population <- NULL
  expect_true(inherits(grouped_gen_tbl, "gen_tbl"))
  expect_false(inherits(grouped_gen_tbl, "grouped_df"))
  expect_false(inherits(grouped_gen_tbl, "grouped_gen_tbl"))

  # removing the genotypes drops the gen_tbl classes
  grouped_gen_tbl <- load_example_gt("grouped_gen_tbl")
  grouped_gen_tbl$genotypes <- NULL
  expect_false(inherits(grouped_gen_tbl, "gen_tbl"))
  expect_false(inherits(grouped_gen_tbl, "grouped_gen_tbl"))
  expect_true(inherits(grouped_gen_tbl, "grouped_df"))
})

test_that("mutate dropping genotypes returns a plain ungrouped tibble", {
  expect_message(
    test_tbl <- test_gt %>% mutate(genotypes = NULL),
    "this is no longer a 'gen_tbl'"
  )
  expect_identical(class(test_tbl), c("tbl_df", "tbl", "data.frame"))
  expect_equal(summarise(test_tbl, n = n())$n, 3)
  expect_message(
    test_tbl <- test_gt %>% mutate(z = 1, .keep = "none"),
    "this is no longer a 'gen_tbl'"
  )
  expect_identical(class(test_tbl), c("tbl_df", "tbl", "data.frame"))
})

test_that("group_by without variables returns an ungrouped gen_tibble", {
  test_ungrouped <- test_gt %>% group_by()
  expect_identical(class(test_ungrouped), class(test_gt))
  expect_equal(count_loci(test_ungrouped), 6)
  expect_identical(
    class(test_ungrouped %>% filter(id != "a")),
    class(test_gt)
  )
  # regrouping a grouped gen_tibble with no variables ungroups it
  test_ungrouped <- test_gt %>%
    group_by(population) %>%
    group_by()
  expect_identical(class(test_ungrouped), class(test_gt))
  gen_tbl_sf <- load_example_gt("gen_tbl_sf")
  expect_identical(class(gen_tbl_sf %>% group_by()), class(gen_tbl_sf))
})

test_that("ungroup keeps the sf and gen_tbl classes", {
  gen_tbl_sf <- load_example_gt("gen_tbl_sf")
  grouped_gen_tbl_sf <- load_example_gt("grouped_gen_tbl_sf")
  ungrouped_sf <- grouped_gen_tbl_sf %>% ungroup()
  expect_identical(class(ungrouped_sf), class(gen_tbl_sf))
  expect_identical(attr(ungrouped_sf, "sf_column"), "geometry")
  expect_true(inherits(sf::st_geometry(ungrouped_sf), "sfc"))
  expect_identical(
    class(grouped_gen_tbl_sf %>% ungroup(population)),
    class(gen_tbl_sf)
  )
  expect_identical(
    class(load_example_gt("grouped_gen_tbl") %>% ungroup()),
    class(test_gt)
  )
  # removing only some of the grouping variables keeps the others
  test_partial <- test_gt %>%
    mutate(k = c(1, 1, 2)) %>%
    group_by(population, k) %>%
    ungroup(k)
  expect_identical(
    class(test_partial),
    c("grouped_gen_tbl", "grouped_df", class(test_gt))
  )
  expect_equal(group_vars(test_partial), "population")
})

test_that("$<- can remove the geometry from a gen_tibble with sf", {
  gen_tbl_sf <- load_example_gt("gen_tbl_sf")
  gen_tbl_sf$geometry <- NULL
  expect_false(inherits(gen_tbl_sf, "sf"))
  expect_false("geometry" %in% names(gen_tbl_sf))
  expect_true(inherits(gen_tbl_sf, "gen_tbl"))

  grouped_gen_tbl_sf <- load_example_gt("grouped_gen_tbl_sf")
  grouped_gen_tbl_sf$geometry <- NULL
  expect_false(inherits(grouped_gen_tbl_sf, "sf"))
  expect_false("geometry" %in% names(grouped_gen_tbl_sf))
  expect_true(inherits(grouped_gen_tbl_sf, "grouped_gen_tbl"))
})
