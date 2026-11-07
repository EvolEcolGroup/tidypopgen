# Combine a gen_tibble to a data.frame or tibble by column

A [`cbind()`](https://rdrr.io/r/base/cbind.html) method to merge
`gen_tibble` objects with data.frames and normal tibbles. Whilst this
works, it is not ideal as it does not check the order of the tables, and
we suggest that you use
[`dplyr::left_join()`](https://dplyr.tidyverse.org/reference/mutate-joins.html)
instead. Note that `cbind` will not combine two `gen_tibbles` (i.e. it
will NOT combine markers for the same individuals)

## Usage

``` r
# S3 method for class 'gen_tbl'
cbind(..., deparse.level = 1)

# S3 method for class 'grouped_gen_tbl'
cbind(..., deparse.level = 1)
```

## Arguments

- ...:

  a gen_tibble and a data.frame or tibble

- deparse.level:

  an integer controlling the construction of column names. See
  [`cbind`](https://rdrr.io/r/base/cbind.html) for details.

## Value

a `gen_tibble`

## Details

As in
[`dplyr::bind_cols()`](https://dplyr.tidyverse.org/reference/bind_cols.html)
(which `dplyr` uses to [`cbind()`](https://rdrr.io/r/base/cbind.html)
grouped tibbles), duplicated column names are made unique (e.g.
`population...2` and `population...4`), with a message. If a grouping
variable is renamed, the result is no longer grouped. The columns `id`
and `genotypes` of the `gen_tibble` can not be duplicated, as renaming
them would break it.

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

# Create a dataframe to combine with the gen_tibble
df <- data.frame(region = c("A", "A", "B", "B", "A", "B", "B"))

# Combine the gen_tibble with the dataframe
example_gt <- cbind(example_gt, df)
```
