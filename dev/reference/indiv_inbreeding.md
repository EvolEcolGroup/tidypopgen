# Individual inbreeding coefficient

This function calculates the inbreeding coefficient for each individual
based on the beta estimate from Weir and Goudet (2017).

## Usage

``` r
indiv_inbreeding(.x, ...)

# S3 method for class 'gen_tbl'
indiv_inbreeding(.x, method = c("WG17"), allele_sharing_mat = NULL, ...)

# S3 method for class 'vctrs_bigSNP'
indiv_inbreeding(.x, method = c("WG17"), allele_sharing_mat = NULL, ...)

# S3 method for class 'grouped_gen_tbl'
indiv_inbreeding(
  .x,
  method = c("WG17"),
  allele_sharing_mat = NULL,
  type = c("vector", "list"),
  ...
)
```

## Arguments

- .x:

  a vector of class `vctrs_bigSNP` (usually the `genotype` column of a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
  object), or a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md).

- ...:

  other arguments passed to specific methods, currently unused.

- method:

  currently only "WG17" (for Weir and Goudet 2017).

- allele_sharing_mat:

  optional and only relevant for "WG17", the matrix of Allele Sharing
  returned by
  [`pairwise_allele_sharing()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pairwise_allele_sharing.md)
  with `as_matrix=TRUE`. As a number of statistics can be derived from
  the Allele Sharing matrix, it is sometimes more efficient to
  pre-compute this matrix. It is not possible to use this with grouped
  tibbles.

- type:

  only for a grouped
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md):
  "vector" (the default), one value per individual, in the order of the
  rows; or "list" (deprecated), one vector per group, in the order of
  [`dplyr::group_keys()`](https://dplyr.tidyverse.org/reference/group_data.html).

## Value

a numeric vector of inbreeding coefficients, one per individual in the
order of the rows, named by the ids of the individuals (when
`allele_sharing_mat` is given, by its row names). For a grouped
[`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md),
the coefficients are computed within each group (with `type = "list"`,
deprecated, they are returned as a list of vectors, one per group; an
empty group gives an empty vector).

## Details

`indiv_inbreeding()` is an S3 generic, with methods for
[`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
objects (grouped or not) and for their `genotypes` column. Other
packages can add methods for their own classes (e.g. tibbles of genotype
likelihoods).

For a grouped
[`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md),
the coefficients are computed within each group (the beta estimate is
relative to the other individuals of the group), giving the same values
as `group_map(.x, ~ indiv_inbreeding(.x))`, but much faster: the
genotypes are read only once, and only the pairs of individuals within
each group are compared. As for an ungrouped
[`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md),
the result has one value per individual, in the order of the rows, so
that it can be added as a column (e.g.
`grouped_gt$inbreeding <- indiv_inbreeding(grouped_gt)`).

The method for the `genotypes` column computes the matrix of allele
sharing of the individuals in the column, so that `indiv_inbreeding()`
can be used within `dplyr` verbs, e.g.
`mutate(inbreeding = indiv_inbreeding(genotypes))`. Within a grouped
`mutate()`, `dplyr` calls it on the individuals of each group in turn,
so it gives the same values as the method for a grouped
[`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md),
but it reads the genotypes once for each group.

## References

Weir, BS and Goudet J (2017) A Unified Characterization of Population
Structure and Relatedness. Genetics (2017) 206:2085

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

example_gt %>% indiv_inbreeding(method = "WG17")
#>           a           b           c           d           e           f 
#> -0.60305344  0.51908397  0.03816794  0.19847328 -0.20229008  0.51908397 
#>           g 
#> -0.44274809 
```
