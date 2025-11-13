# Load example gen_tibble

This function creates a `gen_tibble` object for use in examples in
documentation.

## Usage

``` r
load_example_gt(
  type = c("gen_tbl", "grouped_gen_tbl", "grouped_gen_tbl_sf", "gen_tbl_sf")
)
```

## Arguments

- type:

  a character string indicating the type of `gen_tibble` to create:

  - "gen_tbl": a basic gen_tibble with genotype data and metadata

  - "grouped_gen_tbl": same as "gen_tbl" but grouped by population

  - "grouped_gen_tbl_sf": adds spatial features (longitude/latitude) and
    groups by population

  - "gen_tbl_sf": adds spatial features without grouping

## Value

an example object of the class `gen_tbl`.

## Examples

``` r
# This function creates an example gen_tibble object
example_gt <- load_example_gt("gen_tbl")
```
