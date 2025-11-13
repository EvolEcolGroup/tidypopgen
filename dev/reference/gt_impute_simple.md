# Simple imputation based on allele frequencies

This function provides a very simple imputation algorithm for
`gen_tibble` objects by using the mode, mean or sampling from the allele
frequencies. Each locus is imputed independently (and thus linkage
information is ignored).

## Usage

``` r
gt_impute_simple(x, method = c("mode", "mean0", "random"), n_cores = 1)
```

## Arguments

- x:

  a
  [gen_tibble](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
  with missing data

- method:

  one of

  - 'mode': the most frequent genotype

  - 'mean0': the mean rounded to the nearest integer

  - 'random': randomly sample a genotype based on the observed allele
    frequencies

- n_cores:

  the number of cores to be used

## Value

a
[gen_tibble](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
with imputed genotypes

## Details

This function is a wrapper around
[`bigsnpr::snp_fastImputeSimple()`](https://privefl.github.io/bigsnpr/reference/snp_fastImputeSimple.html).

## See also

[`bigsnpr::snp_fastImputeSimple()`](https://privefl.github.io/bigsnpr/reference/snp_fastImputeSimple.html)
which this function wraps.

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

# Impute the gen_tibble
example_gt <- example_gt %>% gt_impute_simple()

# And we can check it has been imputed
example_gt %>% gt_has_imputed()
#> [1] TRUE
```
