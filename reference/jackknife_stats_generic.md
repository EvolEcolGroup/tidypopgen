# Estimate a block jackknife standard error from per locus estimates

This function is mostly designed for developers: it is a general
function to estimate the genome-wide value of a statistic, together with
its standard error, using a weighted block jackknife. It takes a vector
of per locus estimates, and aggregates them by their mean (or a custom
function), leaving out one block of loci at a time. To compute specific
quantities directly from a `gen_tibble`, use the appropriate
`jackknife_*` functions, e.g
[`jackknife_pairwise_pop_fst()`](https://evolecolgroup.github.io/tidypopgen/reference/jackknife_pairwise_pop_fst.md)
to compute pairwise Fst.

## Usage

``` r
jackknife_stats_generic(
  .x,
  loci_table,
  operator = c("mean", "custom"),
  block_size,
  size_unit = c("snp", "cM", "bp"),
  f = NULL,
  ...
)
```

## Arguments

- .x:

  A vector containing the per locus estimates.

- loci_table:

  a dataframe including at least a column 'chromosome', and additionally
  a column 'position' if `size_unit` is "bp", or a column 'genetic_dist'
  if `size_unit` is "cM".

- operator:

  The operator to use to aggregate the per locus estimates. Either
  "mean" or "custom" to use a custom function `f`.

- block_size:

  The size of the blocks to leave out.

- size_unit:

  Either "snp", "cM" or "bp". If "snp", the block size is in number of
  SNPs. If "cM", the block size is in centiMorgans. If "bp", the block
  size is in base pairs.

- f:

  a custom function to use for the statistic. This function should take
  a vector of locus estimates and return a single value (see Details for
  the functions that are appropriate). Note that, for a custom function,
  the statistic is recomputed from scratch for every block, which is
  slower than using "mean".

- ...:

  Additional arguments to be passed to the custom operator function.

## Value

A tibble with columns: 'estimate' (the statistic computed on all loci),
'jack_estimate' (the bias corrected jackknife estimate), 'se' (the
jackknife standard error), 'z' (the z score, jack_estimate/se), and
'n_blocks' (the number of blocks used).

## Details

The genome is divided into non-overlapping blocks, which never span more
than one chromosome. Blocks can be defined as a number of SNPs, a number
of base pairs, or a size in centiMorgans (which requires the
'genetic_dist' column of the loci table to be filled in). Blocks in base
pairs cover positions 1 to `block_size`, `block_size + 1` to
`2 * block_size`, and so on along each chromosome (as the windows of
[`windows_stats_generic()`](https://evolecolgroup.github.io/tidypopgen/reference/windows_stats_generic.md)
with `step_size` equal to `window_size`), and blocks in centiMorgans
start from a genetic distance of 0. The statistic is then recomputed
leaving out one block at a time, and the standard error is computed with
the weighted delete-m jackknife of Busing et al. (1999), where each
block is weighted by the number of loci with non-missing values that it
contains. Blocks should be larger than the extent of linkage
disequilibrium (e.g. 5 cM, or a few Mb).

The jackknife assumes that the statistic estimates a quantity that does
not change with the number of loci (such as a mean, or a ratio of sums),
and that is a smooth function of the data. A custom function `f` should
thus not be a sum (leaving out a block always makes a sum smaller, which
the jackknife would mistake for a bias), nor a statistic such as the
median or other quantiles, for which the jackknife does not give valid
standard errors. For the sum of the per locus estimates, use
`operator = "mean"`, and multiply the estimates and the standard error
by the number of non-missing values.

## References

Busing FMTA, Meijer E, Van Der Leeden R (1999) Delete-m jackknife for
unequal m. Statistics and Computing 9:3-8.

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

miss_by_locus <- loci_missingness(example_gt)

# Calculate the mean missingness and its standard error,
# using blocks of 2 SNPs
jackknife_stats_generic(miss_by_locus,
  loci_table = show_loci(example_gt),
  operator = "mean", block_size = 2,
  size_unit = "snp"
)
#> # A tibble: 1 × 5
#>   estimate jack_estimate     se     z n_blocks
#>      <dbl>         <dbl>  <dbl> <dbl>    <int>
#> 1   0.0952        0.0952 0.0238     4        3
```
