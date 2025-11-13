# Imputation based XGBoost

This function provides a simple imputation algorithm for `gen_tibble`
objects based on local XGBoost models.

## Usage

``` r
gt_impute_xgboost(
  x,
  alpha = 1e-04,
  size = 200,
  p_train = 0.8,
  n_cor = nrow(x),
  seed = NA,
  n_cores = 1,
  append_error = TRUE
)
```

## Arguments

- x:

  a
  [gen_tibble](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
  with missing data

- alpha:

  Type-I error for testing correlations. Default is `1e-4`.

- size:

  Number of neighbour SNPs to be possibly included in the model imputing
  this particular SNP. Default is `200`.

- p_train:

  Proportion of non missing genotypes that are used for training the
  imputation model while the rest is used to assess the accuracy of this
  imputation model. Default is `0.8`.

- n_cor:

  Number of rows that are used to estimate correlations. Default uses
  them all.

- seed:

  An integer, for reproducibility. Default doesn't use seeds.

- n_cores:

  the number of cores to be used

- append_error:

  boolean, should the xgboost error estimates be appended as an
  attribute to the genotype column of the gen_tibble. If TRUE (the
  default), a matrix of two rows (the number of missing values, and the
  error estimate) and as many columns as the number of loci will be
  appended to the gen_tibble. attr(missing_gt\$genotypes,
  "imputed_errors")

## Value

a
[gen_tibble](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
with imputed genotypes

## Details

This function is a wrapper around
[`bigsnpr::snp_fastImpute()`](https://privefl.github.io/bigsnpr/reference/snp_fastImpute.html).
The error rates from the xgboost, if appended, can be retrieved with
`attr(x$genotypes, "imputed_errors")` where `x` is the `gen_tibble`.

## See also

[`bigsnpr::snp_fastImpute()`](https://privefl.github.io/bigsnpr/reference/snp_fastImpute.html)
which this function wraps.

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

# Impute the gen_tibble
example_gt <- example_gt %>% gt_impute_xgboost()

# And we can check it has been imputed
example_gt %>% gt_has_imputed()
#> [1] TRUE
```
