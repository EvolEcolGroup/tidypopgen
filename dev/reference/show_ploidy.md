# Show the ploidy information of a `gen_tibble`

Extract the ploidy information from a `gen_tibble`. NOTE that this
function does not return the ploidy level for each individual (that is
obtained with
[`indiv_ploidy`](https://evolecolgroup.github.io/tidypopgen/dev/reference/indiv_ploidy.md));
instead, it returns an integer which is either the ploidy level of all
individuals (e.g. 2 indicates all individuals are diploid), or a 0 to
indicate mixed ploidy. The special case of -2 is used to indicate the
presence of pseudo-haploids (i.e. individuals with a ploidy of 2 but for
which we only have information for one allele; the dosages are 0 or 2
for these individuals).

## Usage

``` r
show_ploidy(.x, ...)

# S3 method for class 'tbl_df'
show_ploidy(.x, ...)

# S3 method for class 'vctrs_bigSNP'
show_ploidy(.x, ...)
```

## Arguments

- .x:

  a vector of class `vctrs_bigSNP` (usually the `genotype` column of a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
  object), or a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md).

- ...:

  currently unused.

## Value

the ploidy (0 indicates mixed ploidy)

## See also

[`indiv_ploidy()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/indiv_ploidy.md)

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

example_gt %>% show_ploidy()
#> [1] 2
```
