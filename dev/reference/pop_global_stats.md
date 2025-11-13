# Compute basic population global statistics

This function computes basic population global statistics, following the
notation in Nei 1987 (which in turn is based on Nei and Chesser 1983):

- observed heterozygosity ( \\\hat{h}\_o\\, column header `Ho`)

- expected heterozygosity, also known as gene diversity (
  \\\hat{h}\_s\\, `Hs`)

- total heterozygosity ( \\\hat{h}\_t\\, `Ht`)

- genetic differentiation between subpopulations (\\D\_{st}\\, `Dst`)

- corrected total population diversity (\\h'\_t\\, `Htp`)

- corrected genetic differentiation between subpopulations
  (\\D'\_{st}\\, `Dstp`)

- \\\hat{F}\_{ST}\\ (column header, `Fst`)

- corrected \\\hat{F'}\_{ST}\\ (column header `Fstp`)

- \\\hat{F}\_{IS}\\ (column header, `Fis`)

- Jost's \\\hat{D}\\ (column header, `Dest`)

## Usage

``` r
pop_global_stats(.x, by_locus = FALSE, n_cores = bigstatsr::nb_cores())
```

## Arguments

- .x:

  a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
  (usually grouped, as obtained by using
  [`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html);
  use on a single population will return a number of quantities as
  NA/NaN)

- by_locus:

  boolean, determining whether the statistics should be returned by
  locus(TRUE), or as a single genome wide value (FALSE, the default).

- n_cores:

  number of cores to be used, it defaults to
  [`bigstatsr::nb_cores()`](https://privefl.github.io/bigstatsr/reference/reexports.html)

## Value

a tibble of population statistics, with populations as rows and
statistics as columns

## Details

We use the notation of Nei 1987. That notation was for loci with \\m\\
alleles, but in our case we only have two alleles, so `m=2`.

- Within population observed heterozygosity \\\hat{h}\_o\\ for a locus
  with \\m\\ alleles is defined as:  
  \\\hat{h}\_o= 1-\sum\_{k=1}^{s} \sum\_{i=1}^{m} \hat{X}\_{kii}/s\\  
  where  
  \\\hat{X}\_{kii}\\ represents the proportion of homozygote \\i\\ in
  the sample for the \\k\\th population and  
  \\s\\ the number of populations,  
  following equation 7.38 in Nei(1987) on pp.164.  

- Within population expected heterozygosity (gene diversity)
  \\\hat{h}\_s\\ for a locus with \\m\\ alleles is defined as:  
  \\\hat{h}\_s=(\tilde{n}/(\tilde{n}-1))\[1-\sum\_{i=1}^{m}\bar{\hat{x}\_i^2}-\hat{h}\_o/2\tilde{n}\]\\  
  \#nolint where  
  \\\tilde{n}=s/\sum_k 1/n_k\\ (i.e the harmonic mean of \\n_k\\) and  
  \\\bar{\hat{x}\_i^2}=\sum_k \hat{x}\_{ki}^2/s\\  
  following equation 7.39 in Nei(1987) on pp.164.

- Total heterozygosity (total gene diversity) \\\hat{h}\_t\\ for a locus
  with \\m\\ alleles is defined as:  
  \\\hat{h}\_t = 1-\sum\_{i=1}^{m} \bar{\hat{x}\_i^2} +
  \hat{h}\_s/(\tilde{n}s) - \hat{h}\_o/(2\tilde{n}s)\\  
  where  
  \\\hat{x}\_i=\sum_k \hat{x}\_{ki}/s\\  
  following equation 7.40 in Nei(1987) on pp.164.  

- The amount of gene diversity among samples \\D\_{ST}\\ is defined
  as:  
  \\D\_{ST} = \hat{h}\_t - \hat{h}\_s\\  
  following the equation provided in the text at the top of page 165 in
  Nei(1987).

- The corrected amount of gene diversity among samples \\D'\_{ST}\\ is
  defined as:  
  \\D'\_{ST} = (s/(s-1))D'\_{ST}\\  
  following the equation provided in the text at the top of page 165 in
  Nei(1987).

- Total corrected heterozygosity (total gene diversity) \\\hat{h}\_t\\
  is defined as:  
  \\\hat{h'}\_t = \hat{h}\_s + D'\_{ST}\\  
  following the equation provided in the text at the top of page 165 in
  Nei(1987).

- \\\hat{F}\_{IS}\\ is defined as:  
  \\\hat{F}\_{IS} = 1 - \hat{h}\_o/\hat{h}\_s\\  
  following equation 7.41 in Nei(1987) on pp.164.  

- \\\hat{F}\_{ST}\\ is defined as:  
  \\\hat{F}\_{ST} = 1 - \hat{h}\_s/\hat{h}\_t = D\_{ST}/\hat{h}\_t\\  
  following equation 7.43 in Nei(1987) on pp.165.  

- \\\hat{F'}\_{ST}\\ is defined as:  
  \\\hat{F'}\_{ST} = D'\_{ST}/\hat{h'}\_t\\  
  following the explanation provided in the text at the top of page 165
  in Nei(1987).

- Jost's \\\hat{D}\\ is defined as:  
  \\\hat{D} = (s/(s-1))((\hat{h'}\_t-\hat{h}\_s)/(1-\hat{h}\_s))\\  
  as defined by Jost(2008)

  All these statistics are first computed by locus, and then averaged
  across loci (including any monomorphic locus) to obtain genome-wide
  values. The function uses the same algorithm as
  [`hierfstat::basic.stats()`](https://rdrr.io/pkg/hierfstat/man/basic.stats.html)
  but is optimized for speed and memory usage.

## References

Nei M, Chesser R (1983) Estimation of fixation indexes and gene
diversities. Annals of Human Genetics, 47, 253-259.

Nei M. (1987) Molecular Evolutionary Genetics. Columbia University
Press, pp. 164-165.

Jost L (2008) GST and its relatives do not measure differentiation.
Molecular Ecology, 17, 4015-4026.

## See also

[`hierfstat::basic.stats()`](https://rdrr.io/pkg/hierfstat/man/basic.stats.html)

## Examples

``` r
example_gt <- load_example_gt("grouped_gen_tbl")

# Compute population global statistics
example_gt %>% pop_global_stats()
#>          Ho          Hs          Ht         Dst         Htp        Dstp 
#>  0.44444444  0.40515873  0.42826279  0.02310406  0.43981481  0.03465608 
#>         Fst        Fstp         Fis        Dest 
#>  0.05394832  0.07879699 -0.09696376  0.05826106 

# To return by locus, set by_locus = TRUE
example_gt %>% pop_global_stats(by_locus = TRUE)
#>          Ho        Hs        Ht          Dst       Htp         Dstp         Fst
#> 1 0.3888889 0.3361111 0.5194444  0.183333333 0.6111111  0.275000000  0.35294118
#> 2 0.3888889 0.3968254 0.4656085  0.068783069 0.5000000  0.103174603  0.14772727
#> 3 0.1666667 0.1583333 0.1638889  0.005555556 0.1666667  0.008333333  0.03389831
#> 4 0.5000000 0.3333333 0.3888889  0.055555556 0.4166667  0.083333333  0.14285714
#> 5 0.6111111 0.6984127 0.5661376 -0.132275132 0.5000000 -0.198412698 -0.23364486
#> 6 0.6111111 0.5079365 0.4656085 -0.042328042 0.4444444 -0.063492063 -0.09090909
#>         Fstp         Fis        Dest
#> 1  0.4500000 -0.15702479  0.41422594
#> 2  0.2063492  0.02000000  0.17105263
#> 3  0.0500000 -0.05263158  0.00990099
#> 4  0.2000000 -0.50000000  0.12500000
#> 5 -0.3968254  0.12500000 -0.65789474
#> 6 -0.1428571 -0.20312500 -0.12903226
```
