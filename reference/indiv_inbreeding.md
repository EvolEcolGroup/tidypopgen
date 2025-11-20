# Individual inbreeding coefficient

This function calculates the inbreeding coefficient for each individual
based on the beta estimate from Weir and Goudet (2017).

## Usage

``` r
indiv_inbreeding(.x, method = c("WG17"), allele_sharing_mat = NULL, ...)

# S3 method for class 'tbl_df'
indiv_inbreeding(.x, method = c("WG17"), allele_sharing_mat = NULL, ...)

# S3 method for class 'vctrs_bigSNP'
indiv_inbreeding(.x, method = c("WG17"), allele_sharing_mat = NULL, ...)

# S3 method for class 'grouped_df'
indiv_inbreeding(.x, method = c("WG17"), allele_sharing_mat = NULL, ...)
```

## Arguments

- .x:

  a vector of class `vctrs_bigSNP` (usually the `genotype` column of a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)
  object), or a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md).

- method:

  currently only "WG17" (for Weir and Goudet 2017).

- allele_sharing_mat:

  optional and only relevant for "WG17", the matrix of Allele Sharing
  returned by
  [`pairwise_allele_sharing()`](https://evolecolgroup.github.io/tidypopgen/reference/pairwise_allele_sharing.md)
  with `as_matrix=TRUE`. As a number of statistics can be derived from
  the Allele Sharing matrix, it is sometimes more efficient to
  pre-compute this matrix. It is not possible to use this with grouped
  tibbles.

- ...:

  currently unused.

## Value

a numeric vector of inbreeding coefficients.

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
