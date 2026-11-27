# Restore the class of a tibble from its genotypes column

Some `dplyr` operations return tibbles that have lost the class of a
`gen_tibble`, such as the data frames returned by `pick()` or `across()`
within `summarise()`, `reframe()` or `mutate()`, the result of
`reframe()`, or `as_tibble()`. Their `genotypes` column, however, keeps
its class. The generics of `tidypopgen` therefore have fallback methods
for data frames (and grouped data frames), which use
`restore_genotypes_tbl()` to give back the tibble its class, as given by
`genotypes_tbl_class()` for its `genotypes` column, and then dispatch
again.

## Usage

``` r
genotypes_tbl_class(genotypes, ...)

# Default S3 method
genotypes_tbl_class(genotypes, ...)

# S3 method for class 'vctrs_bigSNP'
genotypes_tbl_class(genotypes, ...)

restore_genotypes_tbl(.x, generic = "this function")
```

## Arguments

- genotypes:

  the `genotypes` column of a tibble

- ...:

  currently unused

- .x:

  a data frame or tibble with a `genotypes` column

- generic:

  the name of the generic for which the class is restored, used in the
  error message if the restored class has no method

## Value

`genotypes_tbl_class()` returns the class of the tibble that holds the
`genotypes` column (`"gen_tbl"` for a `vctrs_bigSNP` column), or `NULL`
for a column of an unknown class.

`restore_genotypes_tbl()` returns `.x` as a tibble with the class given
by `genotypes_tbl_class()`, and its grouped version (e.g.
`"grouped_gen_tbl"`) if `.x` is grouped. A `rowwise_df` is ungrouped. It
stops with an error if `.x` already has that class, as the generic would
otherwise dispatch to the fallback method again.

## Details

Packages that build their own tibbles on `tidypopgen` (e.g. for genotype
likelihoods) add a method to `genotypes_tbl_class()` for the class of
their `genotypes` column, so that their own methods of the `tidypopgen`
generics are used after the class is restored.
