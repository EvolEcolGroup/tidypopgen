# Test if the loci table is ordered

This functions checks that all SNPs in a chromosome are adjacent in the
loci table, and that positions (and, unless
`ignore_genetic_dist = TRUE`, genetic distances) are sorted within
chromosomes. Missing positions (or missing genetic distances, if they
are checked) are an error, as the order of those loci is unknown.

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

  logical, if `TRUE` (the default) only the positions are checked, and
  genetic distances are ignored (they can then be missing).

- ...:

  other arguments passed to specific methods.

## Value

`TRUE` if the loci are ordered, `FALSE` otherwise

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

example_gt %>% is_loci_table_ordered()
#> [1] TRUE
```
