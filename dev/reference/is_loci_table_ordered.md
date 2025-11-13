# Test if the loci table is ordered

This functions checks that all SNPs in a chromosome are adjacent in the
loci table, and that positions are sorted within chromosomes.

## Usage

``` r
is_loci_table_ordered(
  .x,
  error_on_false = FALSE,
  ignore_genetic_dist = TRUE,
  ...
)

# S3 method for class 'tbl_df'
is_loci_table_ordered(
  .x,
  error_on_false = FALSE,
  ignore_genetic_dist = TRUE,
  ...
)

# S3 method for class 'vctrs_bigSNP'
is_loci_table_ordered(
  .x,
  error_on_false = FALSE,
  ignore_genetic_dist = TRUE,
  ...
)
```

## Arguments

- .x:

  a vector of class `vctrs_bigSNP` (usually the `genotype` column of a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
  object), or a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md).

- error_on_false:

  logical, if `TRUE` an error is thrown if the loci are not ordered.

- ignore_genetic_dist:

  logical, if `TRUE` the physical position is not checked.

- ...:

  other arguments passed to specific methods.

## Value

a logical vector defining which loci are transversions

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

example_gt %>% is_loci_table_ordered()
#> [1] TRUE
```
