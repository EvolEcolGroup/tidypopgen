# Return individual ploidy

Returns the ploidy for each individual.

## Usage

``` r
indiv_ploidy(.x, ...)

# S3 method for class 'tbl_df'
indiv_ploidy(.x, ...)

# S3 method for class 'vctrs_bigSNP'
indiv_ploidy(.x, ...)
```

## Arguments

- .x:

  a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md),
  or a vector of class `vctrs_bigSNP` (usually the `genotype` column of
  a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
  object)

- ...:

  currently unused.

## Value

a vector of ploidy, one per individuals in the
[`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

example_gt %>% indiv_ploidy()
#> [1] 2 2 2 2 2 2 2
```
