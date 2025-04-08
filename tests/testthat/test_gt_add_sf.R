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
  test_gt_geom <- test_gt %>% mutate(my_geometry = sf::st_geometry( sf::st_as_sf(
        x = test_indiv_meta,
        coords = c("longitude", "latitude"),
        crs = 4326
      )
    ))
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
}
          )
