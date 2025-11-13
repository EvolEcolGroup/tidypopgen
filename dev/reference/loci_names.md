# Get the names of loci in a `gen_tibble`

Extract the loci names from a `gen_tibble` (or directly from its
`genotype` column).

## Usage

``` r
loci_names(.x, .col = "genotypes", ...)

# S3 method for class 'tbl_df'
loci_names(.x, .col = "genotypes", ...)

# S3 method for class 'vctrs_bigSNP'
loci_names(.x, .col = "genotypes", ...)
```

## Arguments

- .x:

  a vector of class `vctrs_bigSNP` (usually the `genotype` column of a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
  object), or a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md).

- .col:

  the column to be used when a tibble (or grouped tibble is passed
  directly to the function). This defaults to "genotypes" and can only
  take that value. There is no need for the user to set it, but it is
  included to resolve certain tidyselect operations.

- ...:

  currently unused.

## Value

a character vector of names

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")
example_gt %>% loci_names()
#> [1] "rs1" "rs2" "rs3" "rs4" "rs5" "rs6"
```
