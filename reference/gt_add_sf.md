# Add an simple feature geometry to a `gen_tibble`

`gt_add_sf` adds an active sf geometry column to a `gen_tibble` object.
The resulting `gen_tbl` inherits from `sf` and can be used with
functions from the `sf` package. It is possible to either create a
[`sf::sfc`](https://r-spatial.github.io/sf/reference/sfc.html) geometry
column from coordinates, or to provide an existing geometry column
(which will then become the active geometry for `sf`).

## Usage

``` r
gt_add_sf(x, coords = NULL, crs = NULL, sfc_column = NULL)
```

## Arguments

- x:

  a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)
  object

- coords:

  a vector of length 2, giving the names of the x and y columns in `x`
  (i.e. the coordinates, e.g. longitude and latitude). If `coords` is
  not provided, the geometry column must be provided.

- crs:

  the coordinate reference system of the coordinates. If this is not
  set, it will be set to the default value of `sf::st_crs(4326)`.

- sfc_column:

  the name of an
  [`sf::sfc`](https://r-spatial.github.io/sf/reference/sfc.html) column
  to be used as the geometry

## Value

a
[`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)
object with an additional geometry column (and thus belonging also to
`sf` class).

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

# Add some coordinates
example_gt <- example_gt %>% mutate(
  longitude = c(0, 0, 2, 2, 0, 2, 2),
  latitude = c(51, 51, 49, 49, 51, 41, 41)
)

# Convert lat and long to sf:
example_gt <- gt_add_sf(x = example_gt, coords = c("longitude", "latitude"))

# Check class
class(example_gt)
#> [1] "gen_tbl"    "sf"         "tbl_df"     "tbl"        "data.frame"
```
