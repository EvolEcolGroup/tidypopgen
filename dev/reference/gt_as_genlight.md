# Convert a `gen_tibble` to a `genlight` object from `adegenet`

This function converts a `gen_tibble` to a `genlight` object from
`adegenet`

## Usage

``` r
gt_as_genlight(x)
```

## Arguments

- x:

  a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md),
  with population coded as 'population'

## Value

a `genlight` object

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

# Convert to genlight
gt_genlight <- example_gt %>% gt_as_genlight()

# Check object class
class(gt_genlight)
#> [1] "genlight"
#> attr(,"package")
#> [1] "adegenet"
```
