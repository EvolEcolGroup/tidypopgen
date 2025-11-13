# Compute the population expected heterozygosity

This function computes expected population heterozygosity (also referred
to as gene diversity, to avoid the potentially misleading use of the
term "expected" in this context), using the formula of Nei (1987).

## Usage

``` r
pop_het_exp(
  .x,
  by_locus = FALSE,
  include_global = FALSE,
  n_cores = bigstatsr::nb_cores()
)

pop_gene_div(
  .x,
  by_locus = FALSE,
  include_global = FALSE,
  n_cores = bigstatsr::nb_cores()
)
```

## Arguments

- .x:

  a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
  (usually grouped, as obtained by using
  [`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html),
  otherwise the full tibble will be considered as belonging to a single
  population).

- by_locus:

  boolean, determining whether Hs should be returned by locus(TRUE), or
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

Within population expected heterozygosity (gene diversity)
\\\hat{h}\_s\\ for a locus with \\m\\ alleles is defined as:  
\\\hat{h}\_s=\tilde{n}/(\tilde{n}-1)\[1-\sum\_{i}^{m}\bar{\hat{x}\_i^2}-\hat{h}\_o/2\tilde{n}\]\\  
\#nolint

where  
\\\tilde{n}=s/\sum_k 1/n_k\\ (i.e the harmonic mean of \\n_k\\) and  
\\\bar{\hat{x}\_i^2}=\sum_k \hat{x}\_{ki}^2/s\\  
following equation 7.39 in Nei(1987) on pp.164. In our specific case,
there are only two alleles, so \\m=2\\. \\\hat{h}\_s\\ at the genome
level for each population is simply the mean of the locus estimates for
each population.

## References

Nei M. (1987) Molecular Evolutionary Genetics. Columbia University Press

## Examples

``` r
example_gt <- load_example_gt("grouped_gen_tbl")

# Compute expected heterozygosity
example_gt %>% pop_het_exp()
#>      pop1      pop2      pop3 
#> 0.4166667 0.4000000 0.2500000 

# To include the global expected heterozygosity, set include_global = TRUE
example_gt %>% pop_het_exp(include_global = TRUE)
#>      pop1      pop2      pop3    global 
#> 0.4166667 0.4000000 0.2500000 0.4051587 

# To return by locus, set by_locus = TRUE
example_gt %>% pop_het_exp(by_locus = TRUE)
#>           pop1 pop2 pop3
#> [1,] 0.5000000  0.5  0.0
#> [2,] 0.5000000  NaN  0.5
#> [3,] 0.0000000  0.0  0.5
#> [4,] 0.5000000  0.5  0.0
#> [5,] 0.6666667  0.5  NaN
#> [6,] 0.3333333  0.5  NaN
```
