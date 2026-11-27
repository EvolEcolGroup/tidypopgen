# Simulate a toy dataset of unlinked loci

Simulates the genotypes of a few populations at unlinked biallelic loci
under the Balding-Nichols model, and returns them as a `gen_tibble`.
This is a **toy simulator** for quick checks of methods on unlinked loci
(e.g. that a PCA separates populations, or how a statistic behaves with
missing genotypes); it is not meant for inference.

## Usage

``` r
gt_simulate_bn(
  n_indiv = 20,
  n_loci = 1000,
  fst = 0.05,
  missing = 0,
  freq_range = c(0.05, 0.95),
  backingfile = NULL
)
```

## Arguments

- n_indiv:

  the number of individuals per population: one value, or one per
  population.

- n_loci:

  the number of loci.

- fst:

  the divergence (\\F\_{ST}\\) of each population from the ancestral
  population: one value, or one per population. The number of
  populations is the longer of `n_indiv` and `fst`.

- missing:

  the probability that a genotype is missing: one value, or one per
  individual (e.g. to give individuals different missingness).

- freq_range:

  the range of the (uniform) distribution of the ancestral allele
  frequencies.

- backingfile:

  the path, without extension, of the backing files; by default a
  temporary file.

## Value

a `gen_tibble` with columns `id`, `population` and `genotypes`.

## Details

Each locus has an ancestral allele frequency \\p\\ drawn uniformly in
`freq_range`. The frequency in population \\k\\ is drawn from a beta
distribution with mean \\p\\ and variance \\F_k p (1-p)\\ (the
Balding-Nichols model, i.e. independent drift of each population from a
common ancestor, a star-shaped tree), and the genotypes of its
individuals are binomial draws of two alleles at that frequency
(Hardy-Weinberg proportions). Genotypes are then set to missing at
random with probability `missing`.

The model has no linkage, mutation, recombination, demographic history
or tree of populations: for realistic data, simulate with the coalescent
(e.g. msprime, which can be run from R with the slendr package).

All loci are placed on chromosome 1, 1,000 bp apart, but they are
independent. Some loci can be monomorphic in the sample. The genotypes
are generated in memory before being written to the backing file, so
keep the datasets small. Use
[`set.seed()`](https://rdrr.io/r/base/Random.html) beforehand for
reproducible datasets.

## References

Balding D.J. and Nichols R.A. (1995) A method for quantifying
differentiation between populations at multi-allelic loci and its
implications for investigating identity and paternity. Genetica 96:
3-12.

Baumdicker F. et al. (2022) Efficient ancestry and mutation simulation
with msprime 1.0. Genetics 220: iyab229.

Petr M., Haller B.C., Ralph P.L. and Racimo F. (2023) slendr: a
framework for spatio-temporal population genomic simulations on
geographic landscapes. Peer Community Journal 3: e121.

## Examples

``` r
set.seed(1)
# three populations of 20 individuals, 500 loci
sim_gt <- gt_simulate_bn(
  n_indiv = 20, n_loci = 500, fst = c(0.02, 0.05, 0.1)
)
sim_gt
#> # A gen_tibble: 500 loci
#> # A tibble:     60 × 3
#>    id    population  genotypes
#>    <chr> <chr>      <vctr_SNP>
#>  1 ind1  pop1        [1,1,...]
#>  2 ind2  pop1        [0,0,...]
#>  3 ind3  pop1        [1,1,...]
#>  4 ind4  pop1        [0,1,...]
#>  5 ind5  pop1        [0,0,...]
#>  6 ind6  pop1        [1,1,...]
#>  7 ind7  pop1        [1,0,...]
#>  8 ind8  pop1        [0,0,...]
#>  9 ind9  pop1        [0,1,...]
#> 10 ind10 pop1        [1,0,...]
#> # ℹ 50 more rows
# three populations with missing genotypes, increasing from 0 to 90%
# across the individuals of each population
sim_missing <- gt_simulate_bn(
  n_indiv = c(20, 20, 20), n_loci = 500, fst = 0.05,
  missing = rep(seq(0, 0.9, length.out = 20), 3)
)
```
