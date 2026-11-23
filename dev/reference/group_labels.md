# Labels of groups, as tidypopgen names them

This function is designed for developers of packages that add methods to
the generics of tidypopgen for grouped data (e.g. for genotype
likelihoods), so that their outputs label the groups as the methods of
tidypopgen do. A group is labelled by the value of its grouping variable
or, with several grouping variables, by their values pasted with "." (as
the names of the pairs of groups, e.g. `fst_pop1.pop2`).

## Usage

``` r
group_labels(keys)
```

## Arguments

- keys:

  a data frame with one row per group and one column per grouping
  variable, as returned by
  [`dplyr::group_keys()`](https://dplyr.tidyverse.org/reference/group_data.html).

## Value

a character vector with the label of each group.

## See also

[`group_columns()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/group_columns.md)
for the columns that identify the groups in a tidy output.

## Examples

``` r
keys <- data.frame(population = c("a", "a"), sex = c("F", "M"))
group_labels(keys)
#> [1] "a.F" "a.M"
```
