# Checks if a `gen_tibble` has been imputed

This function checks if a dataset has been imputed. Note that having
imputation does not mean that the imputed values are used.

## Usage

``` r
gt_has_imputed(x)
```

## Arguments

- x:

  a `gen_tibble`

## Value

boolean TRUE or FALSE depending on whether the dataset has been imputed

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

# The initial gen_tibble contains no imputed values
example_gt %>% gt_has_imputed()
#> [1] FALSE

# Now impute the gen_tibble
example_gt <- example_gt %>% gt_impute_simple()

# And we can check it has been imputed
example_gt %>% gt_has_imputed()
#> [1] TRUE
```
