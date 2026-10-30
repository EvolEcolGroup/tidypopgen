# Estimate Tajima's D for the whole genome

Note that Tajima's D estimates from data that have been filtered or
ascertained can be difficult to interpret. This function should ideally
be used on sequence data prior to filtering.

## Usage

``` r
pop_tajimas_d(.x, n_cores, block_size, ...)

# S3 method for class 'tbl_df'
pop_tajimas_d(
  .x,
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(nrow(.x), 1),
  ...
)

# S3 method for class 'vctrs_bigSNP'
pop_tajimas_d(
  .x,
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(length(.x), 1),
  ...
)

# S3 method for class 'grouped_df'
pop_tajimas_d(
  .x,
  n_cores = bigstatsr::nb_cores(),
  block_size = bigstatsr::block_size(nrow(.x), 1),
  ...
)
```

## Arguments

- .x:

  a vector of class `vctrs_bigSNP` (usually the `genotypes` column of a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
  object), or a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md).

- n_cores:

  number of cores to be used, it defaults to
  [`bigstatsr::nb_cores()`](https://rdrr.io/pkg/bigparallelr/man/nb_cores.html)

- block_size:

  maximum number of loci read at once.

- ...:

  other arguments passed to specific methods, currently unused.

## Value

A single numeric value (Tajima's D D) for the whole data set, `NA` when
the statistic is not defined. For grouped data a list of Tajima's D D
values (one per group) is returned.

## Details

Nucleotide diversity is estimated at each locus with
[`loci_pi()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/loci_pi.md),
so this function works with diploid, polyploid, mixed-ploidy and
pseudohaploid data.

Missing genotypes are accounted for following Bailey et al. (2025), as
implemented in `pixy`: the sample size \\n_x\\ at each locus \\x\\ is
the number of non-missing allele copies at that locus (i.e. the sum of
the ploidy of the genotyped individuals, with pseudohaploids
contributing one copy each). The two estimators of \\\theta\\ are then:

- \\\hat{\theta}\_\pi = \sum_x \pi_x\\, where \\\pi_x\\ is the
  nucleotide diversity at locus \\x\\, estimated from its \\n_x\\ allele
  copies (loci with fewer than 2 allele copies contribute 0);

- \\\hat{\theta}\_W = \sum_x 1 / a\_{n_x}\\, summed over the segregating
  sites, where \\a\_{n_x} = \sum\_{i=1}^{n_x - 1} 1/i\\ (i.e. each
  segregating site is weighted by its own sample size).

The standard deviation of \\\hat{\theta}\_\pi - \hat{\theta}\_W\\ is
computed with the formula of Tajima (1989), using the number of
segregating sites \\S\\ and the mean sample size across loci with at
least one non-missing allele copy (rounded to the nearest integer, as
the formula is only defined for integer sample sizes). When there are no
missing data, all loci have the same sample size and the statistic is
identical to the classic Tajima's D. Loci with no genotyped individuals
are ignored. Tajima's D is `NA` when there are no segregating sites, or
when the mean sample size is smaller than 2.

## References

Tajima F (1989) Statistical method for testing the neutral mutation
hypothesis by DNA polymorphism. Genetics 123:585-595.

Bailey N, Stevison L, Samuk K (2025) Correcting for bias in estimates of
\\\theta_W\\ and Tajima's D from missing data in next-generation
sequencing. Molecular Ecology Resources, e14104.
[doi:10.1111/1755-0998.14104](https://doi.org/10.1111/1755-0998.14104)

## Examples

``` r
example_gt <- load_example_gt("grouped_gen_tbl")

# Compute Tajima's D
example_gt %>% pop_tajimas_d()
#> [[1]]
#> [1] 0.944084
#> 
#> [[2]]
#> [1] -0.780123
#> 
#> [[3]]
#> [1] -Inf
#> 
```
