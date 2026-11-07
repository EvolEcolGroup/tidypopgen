# Compute the Identity by State Matrix for a `gen_tibble` object

This function computes the IBS matrix.

## Usage

``` r
pairwise_ibs(
  x,
  as_matrix = FALSE,
  type = c("proportion", "adjusted_counts", "raw_counts"),
  block_size = bigstatsr::block_size(nrow(x)),
  n_cores = bigstatsr::nb_cores()
)
```

## Arguments

- x:

  a `gen_tibble` object.

- as_matrix:

  boolean, determining whether the results should be a square
  symmetrical matrix (TRUE), or a tidied tibble (FALSE, the default)

- type:

  one of:

  - "proportion" (the default): the proportion of alleles that are
    identical by state, a similarity (as `--distance ibs flat-missing`
    in PLINK, which writes it to a .mibs file).

  - "adjusted_counts": the allele count distance, a dissimilarity (as
    `--distance flat-missing` in PLINK, which writes it to a .dist
    file): the number of alleles that differ between two individuals,
    scaled up by the number of loci over the number of loci genotyped in
    both individuals, to correct for missing genotypes. It equals
    `2 * n_loci * (1 - proportion)`, so it is 0 for identical
    individuals and at most `2 * n_loci`. PLINK's default correction for
    missing genotypes (without `flat-missing`) is not available (see
    Details).

  - "raw_counts": the counts of identical alleles and of non-missing
    alleles, from which the two other quantities are computed.

- block_size:

  maximum number of loci read at once. Genotypes are compared in a
  compact bit-packed form, so this has little effect on memory use.

- n_cores:

  number of cores to be used, it defaults to
  [`bigstatsr::nb_cores()`](https://rdrr.io/pkg/bigparallelr/man/nb_cores.html)

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
[`gt_set_imputed()`](https://evolecolgroup.github.io/tidypopgen/reference/gt_set_imputed.md)),
and treated as missing otherwise.

Missing genotypes are handled as in PLINK with the `flat-missing`
modifier (`--distance ibs flat-missing` for "proportion",
`--distance flat-missing` for "adjusted_counts"): each pair of
individuals is compared only on the loci genotyped in both, and all loci
have the same weight. By default (without `flat-missing`), PLINK 1.9
instead divides each distance by one minus the sum of the average
contributions to the distance of the loci missing in the pair, which
depend on their allele frequencies. This default correction is not
available in `pairwise_ibs()`, so, when genotypes are missing, its
values are only comparable to those of PLINK run with `flat-missing`.
Without missing genotypes, both corrections give the same values.

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
```
