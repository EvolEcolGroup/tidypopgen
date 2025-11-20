# Sets a `gen_tibble` to use imputed data

This function sets or unsets the use of imputed data. For some analysis,
such as PCA, that does not allow for missing data, we have to use
imputation, but for other analysis it might be preferable to allow for
missing data.

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
