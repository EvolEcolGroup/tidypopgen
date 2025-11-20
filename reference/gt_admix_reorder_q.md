# Reorder the q matrices based on the grouping variable

This function reorders the q matrices in a `gt_admix` object based on
the grouping variable. This is useful before plotting when the samples
from each group are not adjacent to each other in the q matrix.

## Usage

``` r
gt_admix_reorder_q(x, group = NULL)
```

## Arguments

- x:

  a `gt_admix` object, possibly with a grouping variable

- group:

  a character vector with the grouping variable (if there is no grouping
  variable info in `x`)

## Value

a `gt_admix` object with the q matrices reordered

## Examples

``` r
# run the example only if we have the package installed
if (requireNamespace("LEA", quietly = TRUE)) {
  example_gt <- load_example_gt("gen_tbl")

  # Create a gt_admix object
  admix_obj <- example_gt %>% gt_snmf(k = 1:3, project = "force")

  # The $id in admix_obj is the same as in the gen_tibble
  admix_obj$id

  # Reorder the q matrices based on the grouping variable
  admix_obj <- gt_admix_reorder_q(admix_obj,
    group = example_gt$population
  )

  # The $id in admix_obj is now reordered according to the population
  admix_obj$id
}
#> [1] "a" "b" "e" "c" "d" "f" "g"
```
