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
[`bigsnpr::snp_fastImputeSimple()`](https://privefl.github.io/bigsnpr/reference/snp_fastImputeSimple.html),
which imputes the whole file backing matrix (FBM). The allele
frequencies used to impute each locus are therefore estimated from all
the individuals in the FBM, so the `gen_tibble` must contain all of them
(in any order). If some individuals have been removed (e.g. with
[`filter()`](https://rdrr.io/r/stats/filter.html)), use
[`gt_update_backingfile()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_update_backingfile.md)
to resave the FBM with only the individuals in the `gen_tibble` before
imputing.

As each locus is imputed independently, the `gen_tibble` does not need
to contain all the loci in the FBM: the loci that have been removed
(e.g. with
[`select_loci()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/select_loci.md))
are also imputed, but they do not affect the imputation of the loci in
the `gen_tibble`. If only a few of the loci in the FBM are kept,
resaving the FBM with
[`gt_update_backingfile()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_update_backingfile.md)
first avoids the time spent imputing loci that are not used.

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
