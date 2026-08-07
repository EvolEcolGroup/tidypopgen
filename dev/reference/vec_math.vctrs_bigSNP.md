# Mathematical predicates for vctrs_bigSNP vectors

Provides support for the mathematical predicate functions
[`is.finite()`](https://rdrr.io/r/base/is.finite.html),
[`is.infinite()`](https://rdrr.io/r/base/is.finite.html), and
[`is.nan()`](https://rdrr.io/r/base/is.finite.html) on `vctrs_bigSNP`
vectors.

## Usage

``` r
vec_math.vctrs_bigSNP(.fn, .x, ...)
```

## Arguments

- .fn:

  Name of the mathematical function being dispatched.

- .x:

  A `vctrs_bigSNP` vector.

- ...:

  Unused additional arguments.

## Value

A logical vector for
[`is.finite()`](https://rdrr.io/r/base/is.finite.html),
[`is.infinite()`](https://rdrr.io/r/base/is.finite.html), and
[`is.nan()`](https://rdrr.io/r/base/is.finite.html).

## Details

These methods are primarily required for compatibility with RStudio's
object inspector and other tooling that may call mathematical predicates
on custom vctrs classes.

All other Math-group generics are deliberately unsupported and will
raise an error.
