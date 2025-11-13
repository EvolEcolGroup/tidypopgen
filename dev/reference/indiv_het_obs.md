# Estimate individual observed heterozygosity

Estimate observed heterozygosity (H_obs) for each individual (i.e. the
frequency of loci that are heterozygous in an individual).

## Usage

``` r
indiv_het_obs(.x, as_counts = FALSE, ...)

# S3 method for class 'tbl_df'
indiv_het_obs(.x, as_counts = FALSE, ...)

# S3 method for class 'vctrs_bigSNP'
indiv_het_obs(.x, as_counts = FALSE, ...)
```

## Arguments

- .x:

  a vector of class `vctrs_bigSNP` (usually the `genotype` column of a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
  object), or a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md).

- as_counts:

  logical, if `TRUE`, return a matrix with two columns: the number of
  heterozygotes and the number of missing values for each individual.
  These quantities can be useful to compute more complex quantities.

- ...:

  currently unused.

## Value

either:

- a vector of heterozygosities, one per individuals in the
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)

- a matrix with two columns, where the first is the number of
  heterozygous loci for each individual and the second is the number of
  missing values for each individual

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

example_gt %>% indiv_het_obs()
#> [1] 0.6666667 0.2000000 0.4000000 0.3333333 0.5000000 0.2000000 0.6000000

# For observed heterozygosity as counts:
example_gt %>% indiv_het_obs(as_counts = TRUE)
#>      het_n na_n
#> [1,]     4    0
#> [2,]     1    1
#> [3,]     2    1
#> [4,]     2    0
#> [5,]     3    0
#> [6,]     1    1
#> [7,]     3    1
```
