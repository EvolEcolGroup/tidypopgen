# Checks if a `gen_tibble` uses imputed data

This function checks if a dataset uses imputed data. Note that it is
possible to have a dataset that has been imputed but it is currently not
using imputation.

## Usage

``` r
gt_uses_imputed(x)
```

## Arguments

- x:

  a `gen_tibble`

## Value

boolean TRUE or FALSE depending on whether the dataset is using the
imputed values

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

# Impute the gen_tibble
example_gt <- example_gt %>% gt_impute_simple()

# Check whether the gen_tibble uses imputed values
example_gt %>% gt_uses_imputed()
#> [1] FALSE
```
