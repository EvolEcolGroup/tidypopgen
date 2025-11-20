# Count the number of loci in a `gen_tibble`

Count the number of loci in `gen_tibble` (or directly from its
`genotype` column).

## Usage

``` r
count_loci(.x, ...)

# S3 method for class 'tbl_df'
count_loci(.x, ...)

# S3 method for class 'vctrs_bigSNP'
count_loci(.x, ...)
```

## Arguments

- .x:

  a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md),
  or a vector of class `vctrs_bigSNP` (usually the `genotype` column of
  a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)
  object).

- ...:

  currently unused.

## Value

the number of loci

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

example_gt %>% count_loci()
#> [1] 6
```
