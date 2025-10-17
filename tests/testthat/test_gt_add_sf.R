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
  population = c("pop1", "pop1", "pop2", "pop2", "pop1", "pop3", "pop3"),
  longitude = c(0, 0, 2, 2, 0, 2, 2),
  latitude = c(51, 51, 49, 49, 51, 41, 41)
)
test_loci <- data.frame(
  name = paste0("rs", 1:6),
  chromosome = paste0("chr", c(1, 1, 1, 1, 2, 2)),
  position = as.integer(c(3, 5, 65, 343, 23, 456)),
  genetic_dist = as.double(rep(0, 6)),
  allele_ref = c("A", "T", "C", "G", "C", "T"),
  allele_alt = c("T", "C", NA, "C", "G", "A")
)

test_that("gt_add_sf() works with gen_tibble", {
  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  )
  test_gt_from_sf <- gt_add_sf(
    x = test_gt,
    coords = c("longitude", "latitude"),
  )
  # this is both a gen_tibble and an sf object
  expect_true(inherits(test_gt_from_sf, "gen_tbl"))
  expect_true(inherits(test_gt_from_sf, "sf"))
  # now filter by population
  test_gt_from_sf <- test_gt_from_sf %>%
    dplyr::filter(population == "pop1")
  # check that we only have 3 rows
  expect_equal(nrow(test_gt_from_sf), 3)
  # and that the population is pop1
  expect_equal(unique(test_gt_from_sf$population), "pop1")
  # and that we still inherit from gen_tibble
  expect_true(inherits(test_gt_from_sf, "gen_tbl"))
  expect_true(inherits(test_gt_from_sf, "sf"))
  # now drop the population column
  test_gt_from_sf <- test_gt_from_sf %>%
    dplyr::select(-population)
  # check that we inherit from gen_tibble
  expect_true(inherits(test_gt_from_sf, "gen_tbl"))
  expect_true(inherits(test_gt_from_sf, "sf"))

  # now check that the geometry column is sticky
  test_gt_from_sf <- test_gt_from_sf %>%
    dplyr::select(-geometry)
  # check that we still have a geometry column
  expect_true("geometry" %in% names(test_gt_from_sf))
  # check that we inherit from gen_tibble
  expect_true(inherits(test_gt_from_sf, "gen_tbl"))
  expect_true(inherits(test_gt_from_sf, "sf"))
  # now drop the geometry column
  test_gt_minus_sf <- test_gt_from_sf %>%
    sf::st_drop_geometry()
  # check that we don't have a geometry column
  expect_false("geometry" %in% names(test_gt_minus_sf))
  # check that we inherit from gen_tibble
  expect_true(inherits(test_gt_minus_sf, "gen_tbl"))
  expect_false(inherits(test_gt_minus_sf, "sf"))

  # now add the geometry column manually and then set it as active
  test_gt_geom <- test_gt %>% mutate(my_geometry = sf::st_geometry(sf::st_as_sf(
    x = test_indiv_meta,
    coords = c("longitude", "latitude"),
    crs = 4326
  )))
  # we now have a myt_geometry column
  expect_true("my_geometry" %in% names(test_gt_geom))
  # and it is an sfc column
  expect_true(inherits(test_gt_geom$my_geometry, "sfc"))
  # check that we inherit from gen_tibble
  expect_true(inherits(test_gt_geom, "gen_tbl"))
  expect_false(inherits(test_gt_geom, "sf"))
  # now set it as active
  test_gt_geom <- test_gt_geom %>%
    gt_add_sf(sfc_column = "my_geometry")
  # check that we inherit from gen_tibble
  expect_true(inherits(test_gt_geom, "gen_tbl"))
  expect_true(inherits(test_gt_geom, "sf"))
  # check that we have the correct crs
  expect_equal(sf::st_crs(test_gt_geom), sf::st_crs(4326))
})


test_that("gt_add_sf gives the correct errors", {
  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  )
  # check that we error if coordinates is not of length 2
  expect_error(
    gt_add_sf(
      x = test_gt,
      coords = c("longitude", "latitude", "extra")
    ),
    "coords must be a vector of length 2"
  )
  # error if one of the coordinates is not a column
  expect_error(
    gt_add_sf(
      x = test_gt,
      coords = c("longitude", "extra")
    ),
    "coords must be a vector of length 2, with names of columns in x"
  )
  # error if we give both coordinates and sf_column_name
  expect_error(
    gt_add_sf(
      x = test_gt,
      coords = c("longitude", "latitude"),
      sfc_column = "geometry"
    ),
    "You must provide either coords or sfc_column, not both"
  )
  # error if sfc_column does not exist
  expect_error(
    gt_add_sf(
      x = test_gt,
      sfc_column = "extra"
    ),
    "sfc_column 'extra' does not exist in x"
  )
  # at least one of coords or sfc_column must be provided
  expect_error(
    gt_add_sf(
      x = test_gt
    ),
    "You must provide either coords or sfc_column"
  )
})

test_that("retain sf class after imputing and augmenting pca", {
  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  )
  test_gt_from_sf <- gt_add_sf(
    x = test_gt,
    coords = c("longitude", "latitude"),
  )
  # impute the gt and check class
  test_gt_from_sf_impute <-
    gt_impute_simple(test_gt_from_sf, method = "mode", n_cores = 1)
  expect_equal(class(test_gt_from_sf), class(test_gt_from_sf_impute))
  # augment the gt and check class
  test_pca <- test_gt_from_sf_impute %>% gt_pca_randomSVD(k = 3)
  augmented_gt <- augment(x = test_pca, data = test_gt_from_sf_impute)
  expect_equal(class(augmented_gt), class(test_gt_from_sf_impute))

  # DAPC
  clusters <- gt_cluster_pca(test_pca, k = 3, n_pca = 3)
  test_cluster_best <- gt_cluster_pca_best_k(clusters,
    stat = "BIC",
    criterion = "min", quiet = TRUE
  )
  test_dapc <- test_cluster_best %>% gt_dapc()
  augmented_gt_dapc <- augment(x = test_dapc, data = test_gt_from_sf_impute)
  expect_equal(class(augmented_gt_dapc), class(test_gt_from_sf_impute))
})

test_that("retain sf class after being saved and reloaded", {
  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  )
  test_gt_from_sf <- gt_add_sf(
    x = test_gt,
    coords = c("longitude", "latitude"),
  )
  file <- tempfile()
  file_names <- gt_save(test_gt_from_sf, file = file, quiet = TRUE)
  reloaded_gt <- gt_load(file_names[1])
  expect_equal(class(test_gt_from_sf), class(reloaded_gt))
})

test_that("retain sf class after reordering", {
  test_loci <- data.frame(
    name = paste0("rs", 1:6),
    chromosome = paste0("chr", c(1, 1, 2, 1, 2, 2)),
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
  test_gt_from_sf <- gt_add_sf(
    x = test_gt,
    coords = c("longitude", "latitude"),
  )
  reordered <-
    gt_order_loci(test_gt_from_sf, use_current_table = FALSE, quiet = TRUE)
  expect_equal(class(reordered), class(test_gt_from_sf))
})

test_that("select_loci and select_loci_if retain sf class", {
  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  )
  test_gt_from_sf <- gt_add_sf(
    x = test_gt,
    coords = c("longitude", "latitude"),
  )
  # select_loci
  test_gt_subset <- test_gt_from_sf %>% select_loci(c("rs1", "rs2", "rs3"))
  expect_equal(class(test_gt_from_sf), class(test_gt_subset))
  # select_loci_if
  test_gt_subset_chr2 <-
    test_gt_from_sf %>% select_loci_if(loci_chromosomes(genotypes) == "chr2")
  expect_equal(class(test_gt_from_sf), class(test_gt_subset_chr2))
})

test_that("merging two gen_tibbles with sf", {
  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  )
  test_gt_from_sf <- gt_add_sf(
    x = test_gt,
    coords = c("longitude", "latitude"),
  )

  # create a new gt to merge
  test_indiv_meta <- data.frame(
    id = c("A"),
    population = c("pop1"),
    longitude = c(6),
    latitude = c(51)
  )
  test_genotypes <- rbind(c(2, 1, 0, NA, 0, 0))
  test_gt2 <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  )
  # add sf to the new gt
  test_gt_from_sf2 <- gt_add_sf(
    x = test_gt2,
    coords = c("longitude", "latitude"),
  )
  # merge the two
  sf_gt_merged <- rbind(test_gt_from_sf, test_gt_from_sf2, quiet = TRUE)
  # geometry dropped after merging
  expect_false(inherits(sf_gt_merged, "sf"))
  # we can add it back
  expect_warning(sf_gt_merged <- gt_add_sf(
    x = sf_gt_merged,
    coords = c("longitude", "latitude"),
  ), "will be overwritten")
  expect_true(inherits(sf_gt_merged, "sf"))
})

test_that("cbind gen_tibble and extra data with sf", {
  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  )
  test_gt_from_sf <- gt_add_sf(
    x = test_gt,
    coords = c("longitude", "latitude"),
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
  # geometry dropped after merging
  test_combined_gt <- cbind(test_gt, new_df)
  expect_false(inherits(test_combined_gt, "sf"))
  # we can add it back
  sf_gt_merged <- gt_add_sf(
    x = test_combined_gt,
    coords = c("longitude", "latitude"),
  )
  expect_true(inherits(sf_gt_merged, "sf"))
})

test_that("gt_add_sf works with group_by", {
  # Check class if we group, then add sf
  test_gt <- gen_tibble(
    x = test_genotypes,
    loci = test_loci,
    indiv_meta = test_indiv_meta,
    quiet = TRUE
  )
  test_gt <- test_gt %>%
    dplyr::group_by(population)
  test_gt_group_sf <- gt_add_sf(
    x = test_gt,
    coords = c("longitude", "latitude"),
  )
  expect_equal(
    class(test_gt_group_sf),
    c(
      "grouped_gen_tbl", "grouped_df", "gen_tbl",
      "sf", "tbl_df", "tbl", "data.frame"
    )
  )

  # Check class if we add sf, then group
  test_gt <- test_gt %>% ungroup()
  test_gt_sf <- gt_add_sf(
    x = test_gt,
    coords = c("longitude", "latitude"),
  )
  test_gt_sf_group <- test_gt_sf %>% group_by(population)
  expect_equal(
    class(test_gt_sf_group),
    c(
      "grouped_gen_tbl", "grouped_df", "gen_tbl",
      "sf", "tbl_df", "tbl", "data.frame"
    )
  )

  # The two should end up as the same class
  expect_equal(class(test_gt_group_sf), class(test_gt_sf_group))
})
