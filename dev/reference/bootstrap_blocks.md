# Draw the blocks of bootstrap replicates

This function is designed for developers of packages that add bootstrap
methods to the generics of tidypopgen (e.g. for genotype likelihoods),
so that their replicates resample the same blocks as the methods of
tidypopgen. The loci are divided into non-overlapping blocks that never
span chromosomes, defined as in
[`jackknife_stats_generic()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/jackknife_stats_generic.md)
(a number of SNPs, of base pairs, or a size in centiMorgans), and each
replicate draws as many blocks as there are, with replacement.

## Usage

``` r
bootstrap_blocks(
  loci_table,
  n_bootstrap,
  block_size,
  size_unit = c("snp", "cM", "bp")
)
```

## Arguments

- loci_table:

  a dataframe including at least a column 'chromosome', and additionally
  a column 'position' if `size_unit` is "bp", or a column 'genetic_dist'
  if `size_unit` is "cM", with one row per locus, in the order in which
  the loci are used.

- n_bootstrap:

  the number of replicates.

- block_size:

  the size of the blocks, in the unit of `size_unit`.

- size_unit:

  either "snp", "cM" or "bp". If "snp", the block size is in number of
  loci. If "cM", the block size is in centiMorgans. If "bp", the block
  size is in base pairs.

## Value

a list with `block_id` (the block of each locus, an integer from 1 to
the number of blocks, numbered in their order of appearance) and
`weights` (an integer matrix with the number of times each block, in
columns, is drawn in each replicate, in rows).

## Details

The blocks are drawn with the random number generator of R: set its seed
(e.g. with [`set.seed()`](https://rdrr.io/r/base/Random.html)) for
reproducible replicates. The matrix of weights has one value per
replicate and block; to bound its memory, its size is capped at 1e9
values (option `tidypopgen.sfs_max_bootstrap_weights`).

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")
set.seed(1)
blocks <- bootstrap_blocks(show_loci(example_gt),
  n_bootstrap = 3,
  block_size = 2, size_unit = "snp"
)
blocks$block_id
#> [1] 1 1 2 2 3 3
blocks$weights
#>      [,1] [,2] [,3]
#> [1,]    0    1    2
#> [2,]    1    2    0
#> [3,]    0    3    0
```
