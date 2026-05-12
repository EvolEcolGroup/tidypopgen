# Compute intersection counts from a data frame of logical columns

Groups rows by every unique combination of the specified logical columns
and returns a summary table with one row per observed intersection,
sorted by decreasing count.

## Usage

``` r
compute_intersections(df, sets)
```

## Arguments

- df:

  A data frame containing (at least) the columns named in `sets`.

- sets:

  Character vector of column names to treat as set-membership
  indicators. Each column must be coercible to `logical`.

## Value

A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
one column per set (`TRUE`/`FALSE`), a `count` column (number of rows in
that intersection), and an `intersection_id` integer giving the rank
order (1 = largest intersection).
