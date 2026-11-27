# Sets a `gen_tibble` to use imputed data

This function sets or unsets the use of imputed data. For some analysis
that does not allow for missing data (e.g. PCA with
`gt_pca(method = "random")`, or
[`gt_pca_autoSVD()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_pca_autoSVD.md)),
we have to use imputation, but for other analysis it might be preferable
to allow for missing data (e.g.
[`gt_pca()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_pca.md),
[`loci_ld_clump()`](https://evolecolgroup.github.io/tidypopgen/reference/loci_ld_clump.md)
and
[`loci_ld_prune()`](https://evolecolgroup.github.io/tidypopgen/reference/loci_ld_prune.md)
handle missing genotypes, and use imputed ones only if they are set).

## Usage

``` r
gt_set_imputed(x, set = NULL)
```

## Arguments

- x:

  a `gen_tibble`

- set:

  a boolean defining whether imputed data should be used

## Value

the gen_tibble, invisibly

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

# Impute the gen_tibble
example_gt <- example_gt %>% gt_impute_simple()

# Check whether the gen_tibble uses imputed values
example_gt %>% gt_uses_imputed()
#> [1] FALSE

# Set the gen_tibble to use imputed values
example_gt %>% gt_set_imputed(TRUE)

# And check that the gen_tibble uses imputed values again
example_gt %>% gt_uses_imputed()
#> [1] TRUE
```
