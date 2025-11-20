# Convert a `gen_tibble` to a data.frame compatible with `hierfstat`

This function converts a `gen_tibble` to a data.frame formatted to be
used by `hierfstat` functions.

## Usage

``` r
gt_as_hierfstat(x)
```

## Arguments

- x:

  a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md),
  with population coded as 'population'

## Value

a data.frame with a column 'pop' and further column representing the
genotypes (with alleles recoded as 1 and 2)

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

# Convert to hierfstat format
gt_hierfstat <- example_gt %>% gt_as_hierfstat()

# Check object class
class(gt_hierfstat)
#> [1] "data.frame"
```
