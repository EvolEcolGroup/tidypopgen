# Show the genotypes of a `gen_tibble`

Extract the genotypes (as a matrix) from a `gen_tibble`.

## Usage

``` r
show_genotypes(.x, indiv_indices = NULL, loci_indices = NULL, ...)

# S3 method for class 'tbl_df'
show_genotypes(.x, indiv_indices = NULL, loci_indices = NULL, ...)

# S3 method for class 'vctrs_bigSNP'
show_genotypes(.x, indiv_indices = NULL, loci_indices = NULL, ...)
```

## Arguments

- .x:

  a vector of class `vctrs_bigSNP` (usually the `genotype` column of a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)
  object), or a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md).

- indiv_indices:

  indices of individuals

- loci_indices:

  indices of loci

- ...:

  currently unused.

## Value

a matrix of counts of the alternative alleles (see
[`show_loci()`](https://evolecolgroup.github.io/tidypopgen/reference/show_loci.md))
to extract information on the alleles for those loci from a
[`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md).

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

example_gt %>% show_genotypes()
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    1    1    0    1    1    0
#> [2,]    2    1    0   NA    0    0
#> [3,]    2   NA    0    0    1    1
#> [4,]    1    0    0    1    0    0
#> [5,]    1    2    0    1    2    1
#> [6,]    0    0    0    0   NA    1
#> [7,]    0    1    1    0    1   NA
```
