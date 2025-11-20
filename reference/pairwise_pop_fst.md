# Compute pairwise population Fst

This function computes pairwise Fst. The following methods are
implemented:

- 'Hudson': Hudson's formulation, as derived in Bhatia et al (2013) for
  diploids. This is the only method that can also be used with
  pseudohaploid data.

- 'Nei87' : Fst according to Nei (1987) - includes the correction for
  heterozygosity when computing Ht (it uses the same formulation as in
  [`hierfstat::pairwise.neifst()`](https://rdrr.io/pkg/hierfstat/man/pairwise.neifst.html)),

- 'WC84' : Weir and Cockerham (1984), correcting for missing data (it
  uses the same formulation as in
  [`hierfstat::pairwise.WCfst()`](https://rdrr.io/pkg/hierfstat/man/pairwise.WCfst.html)).

## Usage

``` r
pairwise_pop_fst(
  .x,
  type = c("tidy", "pairwise"),
  by_locus = FALSE,
  by_locus_type = c("tidy", "matrix", "list"),
  method = c("Hudson", "Nei87", "WC84"),
  return_num_dem = FALSE,
  n_cores = bigstatsr::nb_cores()
)
```

## Arguments

- .x:

  a grouped
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)
  (as obtained by using
  [`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html))

- type:

  type of object to return One of "tidy" or "pairwise" for a pairwise
  matrix of populations. Default is "tidy".

- by_locus:

  boolean, determining whether Fst should be returned by locus(TRUE), or
  as a single genome wide value obtained by taking the ratio of the mean
  numerator and denominator (FALSE, the default).

- by_locus_type:

  type of object to return. One of "tidy", "matrix" or "list". Default
  is "tidy".

- method:

  one of 'Hudson', 'Nei87', and 'WC84'

- return_num_dem:

  returns a list of numerators and denominators for each locus. This is
  useful for creating windowed estimates of Fst (as we need to compute
  the mean numerator and denominator within each window). Default is
  FALSE.

- n_cores:

  number of cores to be used, it defaults to
  [`bigstatsr::nb_cores()`](https://privefl.github.io/bigstatsr/reference/reexports.html)

## Value

if `type=tidy`, a tibble of genome-wide pairwise Fst values with each
pairwise combination as a row if "by_locus=FALSE", else a list including
the tibble of genome-wide values as well as a matrix with pairwise Fst
by locus with loci as rows and and pairwise combinations as columns. If
`type=pairwise`, a matrix of genome-wide pairwise Fst values is
returned.

## Details

For all formulae, the genome wide estimate is obtained by taking the
ratio of the mean numerators and denominators over all relevant SNPs.

## References

Bhatia G, Patterson N, Sankararaman S, Price AL. (2013) Estimating and
Interpreting FST: The Impact of Rare Variants. Genome Research,
23(9):1514–1521.

Nei, M. (1987) Molecular Evolutionary Genetics. Columbia University
Press

Weir, B. S., & Cockerham, C. C. (1984). Estimating F-statistics for the
analysis of population structure. Evolution, 38(6): 1358–1370.

## See also

[`hierfstat::pairwise.neifst()`](https://rdrr.io/pkg/hierfstat/man/pairwise.neifst.html)

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

# For a basic global pairwise Fst calculation:
example_gt %>%
  group_by(population) %>%
  pairwise_pop_fst(method = "Nei87")
#> # A tibble: 3 × 3
#>   population_1 population_2  value
#>   <chr>        <chr>         <dbl>
#> 1 pop1         pop2         0.0320
#> 2 pop1         pop3         0.143 
#> 3 pop2         pop3         0.0500

# With a pairwise matrix:
example_gt %>%
  group_by(population) %>%
  pairwise_pop_fst(method = "Nei87", type = "pairwise")
#>           pop1      pop2      pop3
#> pop1        NA 0.0320197 0.1428571
#> pop2 0.0320197        NA 0.0500000
#> pop3 0.1428571 0.0500000        NA

# To calculate Fst by locus:
example_gt %>%
  group_by(population) %>%
  pairwise_pop_fst(method = "Hudson", by_locus = TRUE)
#> $Fst_by_locus
#> # A tibble: 18 × 3
#>    loci  stat_name       value
#>    <chr> <chr>           <dbl>
#>  1 rs1   fst_pop1.pop2  -0.24 
#>  2 rs1   fst_pop1.pop3   0.6  
#>  3 rs1   fst_pop2.pop3   0.667
#>  4 rs2   fst_pop1.pop2   0.6  
#>  5 rs2   fst_pop1.pop3   0.114
#>  6 rs2   fst_pop2.pop3   0    
#>  7 rs3   fst_pop1.pop2 NaN    
#>  8 rs3   fst_pop1.pop3   0    
#>  9 rs3   fst_pop2.pop3   0    
#> 10 rs4   fst_pop1.pop2  -0.167
#> 11 rs4   fst_pop1.pop3   0.333
#> 12 rs4   fst_pop2.pop3   0    
#> 13 rs5   fst_pop1.pop2  -0.1  
#> 14 rs5   fst_pop1.pop3  -0.6  
#> 15 rs5   fst_pop2.pop3  -0.5  
#> 16 rs6   fst_pop1.pop2  -0.25 
#> 17 rs6   fst_pop1.pop3  -0.333
#> 18 rs6   fst_pop2.pop3  -0.5  
#> 
#> $Fst
#> # A tibble: 3 × 3
#>   population_1 population_2  value
#>   <chr>        <chr>         <dbl>
#> 1 pop1         pop2         0.0345
#> 2 pop1         pop3         0.0556
#> 3 pop2         pop3         0     
#> 
```
