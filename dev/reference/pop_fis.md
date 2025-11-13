# Compute population specific FIS

This function computes population specific FIS, using either the
approach of Nei 1987 (with an algorithm equivalent to the one used by
[`hierfstat::basic.stats()`](https://rdrr.io/pkg/hierfstat/man/basic.stats.html))
or of Weir and Goudet 2017 (with an algorithm equivalent to the one used
by
[`hierfstat::fis.dosage()`](https://rdrr.io/pkg/hierfstat/man/fs.dosage.html)).

## Usage

``` r
pop_fis(
  .x,
  method = c("Nei87", "WG17"),
  by_locus = FALSE,
  include_global = FALSE,
  allele_sharing_mat = NULL
)
```

## Arguments

- .x:

  a grouped
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
  (as obtained by using
  [`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html))

- method:

  one of "Nei87" (based on Nei 1987, eqn 7.41) or "WG17" (for Weir and
  Goudet 2017) to compute FIS

- by_locus:

  boolean, determining whether FIS should be returned by locus(TRUE), or
  as a single genome wide value (FALSE, the default). Note that this is
  only relevant for "Nei87", as "WG17" always returns a single value.

- include_global:

  boolean determining whether, besides the population specific
  estimates, a global estimate should be appended. Note that this will
  return a vector of n populations plus 1 (the global value), or a
  matrix with n+1 columns if `by_locus=TRUE`.

- allele_sharing_mat:

  optional and only relevant for "WG17", the matrix of Allele Sharing
  returned by
  [`pairwise_allele_sharing()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pairwise_allele_sharing.md)
  with `as_matrix=TRUE`. As a number of statistics can be derived from
  the Allele Sharing matrix, it is sometimes more efficient to
  pre-compute this matrix.

## Value

a vector of population specific fis (plus the global value if
`include_global=TRUE`)

## References

Nei M. (1987) Molecular Evolutionary Genetics. Columbia University Press
Weir, BS and Goudet J (2017) A Unified Characterization of Population
Structure and Relatedness. Genetics (2017) 206:2085

## See also

[`hierfstat::basic.stats()`](https://rdrr.io/pkg/hierfstat/man/basic.stats.html)
[`hierfstat::fis.dosage()`](https://rdrr.io/pkg/hierfstat/man/fs.dosage.html)

## Examples

``` r
example_gt <- load_example_gt("grouped_gen_tbl")

# Compute FIS using Nei87
example_gt %>% pop_fis(method = "Nei87")
#>       pop1       pop2       pop3 
#> -0.2333333  0.0000000  0.0000000 

# Compute FIS using WG17
example_gt %>% pop_fis(method = "WG17")
#>        pop1        pop2        pop3 
#> -0.12328767  0.08333333 -0.60000000 

# To include the global FIS, set include_global = TRUE
example_gt %>% pop_fis(method = "Nei87", include_global = TRUE)
#>        pop1        pop2        pop3      global 
#> -0.23333333  0.00000000  0.00000000 -0.09696376 

# To return FIS by locus, set by_locus = TRUE
example_gt %>% pop_fis(method = "Nei87", by_locus = TRUE)
#>               pop1 pop2 pop3
#> [1,] -3.333333e-01    0  NaN
#> [2,] -3.333333e-01  NaN    0
#> [3,]           NaN  NaN    0
#> [4,] -1.000000e+00    0  NaN
#> [5,]  5.000000e-01    0  NaN
#> [6,] -4.440892e-16    0  NaN

# To calculate from a pre-computed allele sharing matrix:
allele_sharing_mat <- pairwise_allele_sharing(example_gt, as_matrix = TRUE)
example_gt %>% pop_fis(
  method = "WG17",
  allele_sharing_mat = allele_sharing_mat
)
#>        pop1        pop2        pop3 
#> -0.12328767  0.08333333 -0.60000000 
```
