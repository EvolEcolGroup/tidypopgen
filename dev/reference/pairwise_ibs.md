# Compute the Identity by State Matrix for a `gen_tibble` object

This function computes the IBS matrix.

## Usage

``` r
pairwise_ibs(x, ...)

# S3 method for class 'tbl_df'
pairwise_ibs(
  x,
  as_matrix = FALSE,
  type = c("proportion", "adjusted_counts", "raw_counts"),
  block_size = bigstatsr::block_size(nrow(x)),
  n_cores = bigstatsr::nb_cores(),
  missing_correction = c("flat", "allele_freq"),
  allele_freqs = NULL,
  ...
)

# S3 method for class 'vctrs_bigSNP'
pairwise_ibs(
  x,
  as_matrix = FALSE,
  type = c("proportion", "adjusted_counts", "raw_counts"),
  block_size = bigstatsr::block_size(length(x)),
  n_cores = bigstatsr::nb_cores(),
  missing_correction = c("flat", "allele_freq"),
  allele_freqs = NULL,
  ...
)
```

## Arguments

- x:

  a `gen_tibble` object, or a vector of class `vctrs_bigSNP` (the
  `genotypes` column of a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md),
  whose individuals are then labelled with the names of the vector).

- ...:

  other arguments passed to specific methods, currently unused.

- as_matrix:

  boolean, determining whether the results should be a square
  symmetrical matrix (TRUE), or a tidied tibble (FALSE, the default)

- type:

  one of:

  - "proportion" (the default): the proportion of alleles that are
    identical by state, a similarity (as `--distance ibs flat-missing`
    in PLINK, which writes it to a .mibs file; with
    `missing_correction = "allele_freq"`, as `--distance ibs`, see
    Details).

  - "adjusted_counts": the allele count distance, a dissimilarity (as
    `--distance flat-missing` in PLINK, which writes it to a .dist
    file): the number of alleles that differ between two individuals,
    scaled up by the number of loci over the number of loci genotyped in
    both individuals, to correct for missing genotypes. It equals
    `2 * n_loci * (1 - proportion)`, so it is 0 for identical
    individuals and, with the default `missing_correction = "flat"`, at
    most `2 * n_loci`. With `missing_correction = "allele_freq"`, as
    `--distance` in PLINK (without `flat-missing`), it can exceed
    `2 * n_loci`, see Details.

  - "raw_counts": the counts of identical alleles and of non-missing
    alleles, from which the two other quantities are computed.

- block_size:

  maximum number of loci read at once. Genotypes are compared in a
  compact bit-packed form, so this has little effect on memory use.

- n_cores:

  number of cores to be used, it defaults to
  [`bigstatsr::nb_cores()`](https://rdrr.io/pkg/bigparallelr/man/nb_cores.html)

- missing_correction:

  the correction for missing genotypes, for "proportion" and
  "adjusted_counts": "flat" (the default, as PLINK 1.9 with the
  `flat-missing` modifier) or "allele_freq" (as PLINK 1.9's default,
  which weights loci by their allele frequencies); see Details. It must
  be "flat" for "raw_counts".

- allele_freqs:

  only used if `missing_correction = "allele_freq"`: the frequency of
  the alternate allele at each locus (as PLINK's `--read-freq`). If
  `NULL` (the default), it is computed from all the individuals of `x`
  (using imputed genotypes if `x` is set to use them).

## Value

for "proportion" and "adjusted_counts", a tidy tibble with one row per
pair of individuals (columns `item1`, `item2` and `value`), or a square
matrix if `as_matrix = TRUE`; for "raw_counts", a list of two
[bigstatsr::FBM](https://privefl.github.io/bigstatsr/reference/FBM-class.html)
matrices, `ibs` (the counts of alleles identical by state) and `valid_n`
(the counts of non-missing alleles, i.e. twice the number of loci
genotyped in both individuals).

## Details

Note that monomorphic sites are currently considered (as in PLINK).
Remove monomorphic sites before running `pairwise_ibs()` if this is a
concern.

Imputed genotypes are used if `x` is set to use them (see
[`gt_set_imputed()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_set_imputed.md)),
and treated as missing otherwise.

Each pair of individuals is compared only on the loci genotyped in both.
Two corrections for missing genotypes are available
(`missing_correction`), which give the same values when no genotypes are
missing:

- "flat" (the default), as PLINK 1.9 with the `flat-missing` modifier
  (`--distance ibs flat-missing` for "proportion",
  `--distance flat-missing` for "adjusted_counts"): all loci have the
  same weight, so the number of alleles that differ between individuals
  `i` and `j`, `d_ij`, is scaled up by `n_loci / n_ij`, where `n_ij` is
  the number of loci genotyped in both ("adjusted_counts"), and
  "proportion" is `1 - d_ij / (2 * n_ij)`.

- "allele_freq", as PLINK 1.9's default (`--distance ibs` and
  `--distance`): loci are weighted by their expected contribution to the
  distance, which depends on their allele frequency `p`. Under
  Hardy-Weinberg equilibrium, the expected number of alleles that differ
  between two individuals at a locus is
  `4 * p * (1 - p) * (p^2 - p + 1)`, and PLINK uses
  `w = p * (1 - p) * (p^2 - p + 1)` as the weight of the locus (0 for
  monomorphic loci). "adjusted_counts" is then `d_ij * W / W_ij`, where
  `W` is the sum of the weights of all loci and `W_ij` that of the loci
  genotyped in both individuals, and "proportion" is
  `1 - adjusted_counts / (2 * n_loci)`. As `W / W_ij` depends on which
  loci are genotyped in both individuals, "adjusted_counts" can exceed
  `2 * n_loci` (and "proportion" can be negative) for pairs genotyped in
  both mostly at loci with small weights (loci with a rare allele), as
  in PLINK. For example, with two loci with frequencies 0.1 and 0.5, a
  pair of opposite homozygotes at the first locus, missing at the
  second, has "adjusted_counts" of about 6.58 (more than
  `2 * n_loci = 4`). As in PLINK, the weights are scaled to sum to just
  under `2^32` and rounded to integers (so values are identical to
  PLINK's), `W_ij` excludes each locus missing in either individual or
  in both (its weight is removed once), and monomorphic loci have no
  weight (a missing genotype at a monomorphic locus does not change the
  distance, but the locus still counts in `n_loci`). PLINK 1.9 computes
  the allele frequencies from the founders only (unless `--nonfounders`
  is used), or reads them with `--read-freq`; `pairwise_ibs()` computes
  them from all individuals of `x` (the same as PLINK when all
  individuals are founders, as in data without a pedigree), unless they
  are given with `allele_freqs` (e.g. from the founders, with
  [`loci_alt_freq()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/loci_alt_freq.md)).
  Loci without any genotype have a frequency of 0.5, as in PLINK.
  PLINK's `1-ibs` modifier (`--distance 1-ibs`, written to a `.mdist`
  file) gives `1 - proportion` with either correction.

Pairs of individuals without any locus genotyped in both (or, with
"allele_freq", without any polymorphic locus genotyped in both) have
`NaN` (or `Inf`) values.

`pairwise_ibs()` is an S3 generic. The method for a
[`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)
computes the matrix with the method for its `genotypes` column (of class
`vctrs_bigSNP`), and labels it with the `id` of the individuals;
grouping is ignored. Other packages can add methods for their own
classes of genotypes (e.g. genotype likelihoods).

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

pairwise_ibs(example_gt, type = "proportion")
#> # A tibble: 21 × 3
#>    item1 item2 value
#>    <chr> <chr> <dbl>
#>  1 a     b     0.8  
#>  2 a     c     0.7  
#>  3 a     d     0.833
#>  4 a     e     0.75 
#>  5 a     f     0.6  
#>  6 a     g     0.7  
#>  7 b     c     0.75 
#>  8 b     d     0.8  
#>  9 b     e     0.5  
#> 10 b     f     0.5  
#> # ℹ 11 more rows

# Alternatively, return a matrix
pairwise_ibs(example_gt, type = "proportion", as_matrix = TRUE)
#>           a    b     c         d         e    f     g
#> a 1.0000000 0.80 0.700 0.8333333 0.7500000 0.60 0.700
#> b 0.8000000 1.00 0.750 0.8000000 0.5000000 0.50 0.500
#> c 0.7000000 0.75 1.000 0.6000000 0.7000000 0.75 0.625
#> d 0.8333333 0.80 0.600 1.0000000 0.5833333 0.70 0.500
#> e 0.7500000 0.50 0.700 0.5833333 1.0000000 0.60 0.500
#> f 0.6000000 0.50 0.750 0.7000000 0.6000000 1.00 0.750
#> g 0.7000000 0.50 0.625 0.5000000 0.5000000 0.75 1.000
#> attr(,"class")
#> [1] "pairwise_matrix" "matrix"          "array"          

# Adjust block_size
pairwise_ibs(example_gt, block_size = 2)
#> # A tibble: 21 × 3
#>    item1 item2 value
#>    <chr> <chr> <dbl>
#>  1 a     b     0.8  
#>  2 a     c     0.7  
#>  3 a     d     0.833
#>  4 a     e     0.75 
#>  5 a     f     0.6  
#>  6 a     g     0.7  
#>  7 b     c     0.75 
#>  8 b     d     0.8  
#>  9 b     e     0.5  
#> 10 b     f     0.5  
#> # ℹ 11 more rows

# Change type
pairwise_ibs(example_gt, type = "adjusted_counts")
#> # A tibble: 21 × 3
#>    item1 item2 value
#>    <chr> <chr> <dbl>
#>  1 a     b       2.4
#>  2 a     c       3.6
#>  3 a     d       2  
#>  4 a     e       3  
#>  5 a     f       4.8
#>  6 a     g       3.6
#>  7 b     c       3  
#>  8 b     d       2.4
#>  9 b     e       6  
#> 10 b     f       6  
#> # ℹ 11 more rows
pairwise_ibs(example_gt, type = "raw_counts")
#> $ibs
#> A Filebacked Big Matrix of type 'double' with 7 rows and 7 columns.
#> 
#> $valid_n
#> A Filebacked Big Matrix of type 'double' with 7 rows and 7 columns.
#> 

# PLINK 1.9's default correction for missing genotypes
pairwise_ibs(example_gt,
  type = "adjusted_counts",
  missing_correction = "allele_freq"
)
#> # A tibble: 21 × 3
#>    item1 item2 value
#>    <chr> <chr> <dbl>
#>  1 a     b      2.40
#>  2 a     c      3.75
#>  3 a     d      2   
#>  4 a     e      3   
#>  5 a     f      5.00
#>  6 a     g      3.59
#>  7 b     c      3.15
#>  8 b     d      2.40
#>  9 b     e      5.99
#> 10 b     f      6.30
#> # ℹ 11 more rows
```
