# Compute population specific Fst

This function computes population specific Fst, using the approach in
Weir and Goudet 2017 (as computed by
[`hierfstat::fst.dosage()`](https://rdrr.io/pkg/hierfstat/man/fs.dosage.html)).

## Usage

``` r
pop_fst(.x, include_global = FALSE, allele_sharing_mat = NULL)
```

## Arguments

- .x:

  a grouped
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)
  (as obtained by using
  [`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html))

- include_global:

  boolean determining whether, besides the population specific Fst, a
  global Fst should be appended. Note that this will return a vector of
  n populations plus 1 (the global value)

- allele_sharing_mat:

  optional, the matrix of Allele Sharing returned by
  [`pairwise_allele_sharing()`](https://evolecolgroup.github.io/tidypopgen/reference/pairwise_allele_sharing.md)
  with `as_matrix=TRUE`. As a number of statistics can be derived from
  the Allele Sharing matrix,

## Value

a vector of population specific Fst (plus the global value if
`include_global=TRUE`)

## References

Weir, BS and Goudet J (2017) A Unified Characterization of Population
Structure and Relatedness. Genetics (2017) 206:2085

## See also

[`hierfstat::fst.dosage()`](https://rdrr.io/pkg/hierfstat/man/fs.dosage.html)

## Examples

``` r
example_gt <- load_example_gt("grouped_gen_tbl")

# Compute FIS using Nei87
example_gt %>% pop_fst()
#>       pop1       pop2       pop3 
#> 0.05246079 0.06544078 0.41590049 

# To include the global Fst, set include_global = TRUE
example_gt %>% pop_fst(include_global = TRUE)
#>       pop1       pop2       pop3     global 
#> 0.05246079 0.06544078 0.41590049 0.17793402 

# To calculate from a pre-computed allele sharing matrix:
allele_sharing_mat <- pairwise_allele_sharing(example_gt, as_matrix = TRUE)
example_gt %>% pop_fst(allele_sharing_mat = allele_sharing_mat)
#>       pop1       pop2       pop3 
#> 0.05246079 0.06544078 0.41590049 
```
