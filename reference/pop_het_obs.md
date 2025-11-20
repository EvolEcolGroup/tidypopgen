# Compute the population observed heterozygosity

This function computes population heterozygosity, using the formula of
Nei (1987).

## Usage

``` r
pop_het_obs(
  .x,
  by_locus = FALSE,
  include_global = FALSE,
  n_cores = bigstatsr::nb_cores()
)
```

## Arguments

- .x:

  a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)
  (usually grouped, as obtained by using
  [`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html),
  otherwise the full tibble will be considered as belonging to a single
  population).

- by_locus:

  boolean, determining whether Ho should be returned by locus(TRUE), or
  as a single genome wide value (FALSE, the default).

- include_global:

  boolean determining whether, besides the population specific
  estimates, a global estimate should be appended. Note that this will
  return a vector of n populations plus 1 (the global value), or a
  matrix with n+1 columns if `by_locus=TRUE`.

- n_cores:

  number of cores to be used, it defaults to
  [`bigstatsr::nb_cores()`](https://privefl.github.io/bigstatsr/reference/reexports.html)

## Value

a vector of mean population observed heterozygosities (if
`by_locus=FALSE`), or a matrix of estimates by locus (rows are loci,
columns are populations, `by_locus=TRUE`)

## Details

Within population observed heterozygosity \\\hat{h}\_o\\ for a locus
with \\m\\ alleles is defined as:  
\\\hat{h}\_o= 1-\sum\_{k=1}^{s} \sum\_{i=1}^{m} \hat{X}\_{kii}/s\\  
where  
\\\hat{X}\_{kii}\\ represents the proportion of homozygote \\i\\ in the
sample for the \\k\\th population and  
\\s\\ the number of populations,  
following equation 7.38 in Nei(1987) on pp.164. In our specific case,
there are only two alleles, so \\m=2\\. For population specific
estimates, the sum is done over a single value of \\k\\. \\\hat{h}\_o\\
at the genome level is simply the mean of the locus estimates.

## References

Nei M. (1987) Molecular Evolutionary Genetics. Columbia University Press

## Examples

``` r
example_gt <- load_example_gt("grouped_gen_tbl")

# Compute expected heterozygosity
example_gt %>% pop_het_obs()
#>      pop1      pop2      pop3 
#> 0.5000000 0.3333333 0.5000000 

# To include the global expected heterozygosity, set include_global = TRUE
example_gt %>% pop_het_obs(include_global = TRUE)
#>      pop1      pop2      pop3    global 
#> 0.5000000 0.3333333 0.5000000 0.4444444 

# To return by locus, set by_locus = TRUE
example_gt %>% pop_het_obs(by_locus = TRUE)
#>           pop1 pop2 pop3
#> [1,] 0.6666667  0.5  0.0
#> [2,] 0.6666667  0.0  0.5
#> [3,] 0.0000000  0.0  0.5
#> [4,] 1.0000000  0.5  0.0
#> [5,] 0.3333333  0.5  1.0
#> [6,] 0.3333333  0.5  1.0
```
