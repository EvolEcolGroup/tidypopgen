# Convert a `gen_tibble` to a `genind` object from `adegenet`

This function converts a `gen_tibble` to a `genind` object from
`adegenet`

## Usage

``` r
gt_as_genind(x)
```

## Arguments

- x:

  a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md),
  with population coded as 'population'

## Value

a `genind` object

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

# Convert to genind
gt_genind <- example_gt %>% gt_as_genind()

# Check object class
class(gt_genind)
#> [1] "genind"
#> attr(,"package")
#> [1] "adegenet"
```
