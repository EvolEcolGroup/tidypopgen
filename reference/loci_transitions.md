# Find transitions

Use the loci table to define which loci are transitions

## Usage

``` r
loci_transitions(.x, .col = "genotypes", ...)

# S3 method for class 'tbl_df'
loci_transitions(.x, .col = "genotypes", ...)

# S3 method for class 'vctrs_bigSNP'
loci_transitions(.x, .col = "genotypes", ...)
```

## Arguments

- .x:

  a vector of class `vctrs_bigSNP` (usually the `genotype` column of a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)
  object), or a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md).

- .col:

  the column to be used when a tibble (or grouped tibble is passed
  directly to the function). This defaults to "genotypes" and can only
  take that value. There is no need for the user to set it, but it is
  included to resolve certain tidyselect operations.

- ...:

  other arguments passed to specific methods.

## Value

a logical vector defining which loci are transitions

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")
example_gt %>% loci_transitions()
#> [1] FALSE  TRUE    NA FALSE FALSE FALSE
```
