# Columns that identify groups in a tidy output, as tidypopgen names them

This function is designed for developers of packages that add methods to
the generics of tidypopgen for grouped data (e.g. for genotype
likelihoods), so that their tidy outputs identify the groups as the
methods of tidypopgen do: a single `group` column when there is one
grouping variable, and one column per grouping variable when there are
several (as
[`dplyr::summarise()`](https://dplyr.tidyverse.org/reference/summarise.html)),
keeping their types (e.g. a factor keeps its levels). For outputs on
pairs of groups, each group of the pair has one column per grouping
variable, named with a suffix (e.g. `population_1` and `population_2`,
with `suffix = "_1"` and `"_2"`), also when there is a single grouping
variable.

## Usage

``` r
group_columns(keys, idx, suffix = NULL, other_names = character())
```

## Arguments

- keys:

  a data frame with one row per group and one column per grouping
  variable, as returned by
  [`dplyr::group_keys()`](https://dplyr.tidyverse.org/reference/group_data.html).

- idx:

  the index of the group of each row of the output.

- suffix:

  `NULL` for outputs with a row per group, or a suffix added to the
  names of the grouping variables, for the columns of one group of a
  pair.

- other_names:

  the names of the other columns of the output: the group columns must
  not have any of these names (an error asks to rename the grouping
  variables).

## Value

a tibble with one row per element of `idx`.

## See also

[`group_labels()`](https://evolecolgroup.github.io/tidypopgen/reference/group_labels.md)
for the labels of the groups.

## Examples

``` r
keys <- data.frame(population = c("a", "a", "b"), sex = c("F", "M", "F"))
group_columns(keys, c(1, 3))
#> # A tibble: 2 × 2
#>   population sex  
#>   <chr>      <chr>
#> 1 a          F    
#> 2 b          F    
group_columns(keys[1], c(1, 3))
#> # A tibble: 2 × 1
#>   group
#>   <chr>
#> 1 a    
#> 2 b    
group_columns(keys, c(1, 3), suffix = "_1")
#> # A tibble: 2 × 2
#>   population_1 sex_1
#>   <chr>        <chr>
#> 1 a            F    
#> 2 b            F    
```
