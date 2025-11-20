# Estimate individual missingness

Estimate missingness for each individual (i.e. the frequency of missing
genotypes in an individual).

## Usage

``` r
indiv_missingness(.x, as_counts, block_size, ...)

# S3 method for class 'tbl_df'
indiv_missingness(
  .x,
  as_counts = FALSE,
  block_size = bigstatsr::block_size(nrow(.x), 1),
  ...
)

# S3 method for class 'vctrs_bigSNP'
indiv_missingness(
  .x,
  as_counts = FALSE,
  block_size = bigstatsr::block_size(length(.x), 1),
  ...
)
```

## Arguments

- .x:

  a vector of class `vctrs_bigSNP` (usually the `genotype` column of a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)
  object), or a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md).

- as_counts:

  boolean defining whether the count of NAs (rather than the rate)
  should be returned. It defaults to FALSE (i.e. rates are returned by
  default).

- block_size:

  maximum number of loci read at once.

- ...:

  currently unused.

## Value

a vector of missingness, one per individuals in the
[`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

example_gt %>% indiv_missingness()
#> [1] 0.0000000 0.1666667 0.1666667 0.0000000 0.0000000 0.1666667 0.1666667

# For missingness as counts:
example_gt %>% indiv_missingness(as_counts = TRUE)
#> [1] 0 1 1 0 0 1 1
```
