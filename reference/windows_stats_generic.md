# Estimate window statistics from per locus estimates

This function is mostly designed for developers: it is a general
function to estimate window statistics from per locus estimates. This
function takes a vector of per locus estimates, and aggregates them by
sum or mean per window. To compute specific quantities directly from a
`gen_tibble`, use the appropriate `window_*` functions, e.g
[`windows_pairwise_pop_fst()`](https://evolecolgroup.github.io/tidypopgen/reference/windows_pairwise_pop_fst.md)
to compute pairwise Fst.

## Usage

``` r
windows_stats_generic(
  .x,
  loci_table,
  operator = c("mean", "sum", "custom"),
  window_size,
  step_size,
  size_unit = c("snp", "bp"),
  min_loci = 1,
  complete = FALSE,
  f = NULL,
  ...
)
```

## Arguments

- .x:

  A vector containing the per locus estimates.

- loci_table:

  a dataframe including at least a column 'chromosome', and additionally
  a column 'position' if `size_unit` is "bp".

- operator:

  The operator to use for the window statistics. Either "mean", "sum" or
  "custom" to use a custom function `.f`.

- window_size:

  The size of the window to use for the estimates.

- step_size:

  The step size to use for the windows.

- size_unit:

  Either "snp" or "bp". If "snp", the window size and step size are in
  number of SNPs. If "bp", the window size and step size are in base
  pairs.

- min_loci:

  The minimum number of loci required to calculate a window statistic.
  If the number of loci in a window is less than this, the window
  statistic will be NA.

- complete:

  Should the function be evaluated on complete windows only? If FALSE,
  the default, then partial computations will be allowed at the end of
  the chromosome.

- f:

  a custom function to use for the window statistics. This function
  should take a vector of locus estimates and return a single value.

- ...:

  Additional arguments to be passed to the custom operator function.

## Value

A tibble with columns: 'chromosome', 'start', 'end', 'stats', and
'n_loci'. The 'stats' column contains the mean of the per locus
estimates in the window, and 'n_loci' contains the number of loci in the
window.

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

miss_by_locus <- loci_missingness(example_gt)

# Calculate mean missingness across windows
windows_stats_generic(miss_by_locus,
  loci_table = show_loci(example_gt),
  operator = "mean", window_size = 1000,
  step_size = 1000, size_unit = "bp",
  min_loci = 1, complete = FALSE
)
#> # A tibble: 2 × 5
#>   chromosome start   end   stat n_loci
#>   <chr>      <dbl> <dbl>  <dbl>  <dbl>
#> 1 chr1           1  1000 0.0714      4
#> 2 chr2           1  1000 0.143       2
```
