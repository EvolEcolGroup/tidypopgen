# Compute Tajima's D for a sliding window

This function computes Tajima's D for a sliding window across each
chromosome.

## Usage

``` r
windows_pop_tajimas_d(
  .x,
  type = c("matrix", "tidy", "list"),
  window_size,
  step_size,
  size_unit = c("snp", "bp"),
  min_loci = 1,
  complete = FALSE
)
```

## Arguments

- .x:

  a (potentially grouped) `gen_tibble` object

- type:

  type of object to return, if using grouped method. One of "matrix",
  "tidy", or "list". Default is "matrix".

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

## Value

if data is not grouped, a data frame with the following columns:

- `chromosome`: the chromosome for the window

- `start`: the starting locus of the window

- `end`: the ending locus of the window

- `tajimas_d`: the Tajima's D for the population if data are grouped,
  either: a data frame as above with the following columns:

- `chromosome`: the chromosome for the window

- `start`: the starting locus of the window

- `end`: the ending locus of the window

- `n_loci`: the number of loci in the window

- `group`: the Tajima's D for the group for the given window (there will
  be as many of these columns as groups in the gen_tibble, and they will
  be named by the grouping levels) a tidy tibble with the following
  columns:

- `chromosome`: the chromosome for the window

- `start`: the starting locus of the window

- `end`: the ending locus of the window

- `n_loci`: the number of loci in the window

- `group`: the name of the group

- `stat`: the Tajima's D for the given group at the given window or a
  list of data frames, one per group, with the following columns:

- `chromosome`: the chromosome for the window

- `start`: the starting locus of the window

- `end`: the ending locus of the window

- `stat`: the Tajima's D for the given window

- `n_loci`: the number of loci in the window

## Examples

``` r
example_gt <- load_example_gt("grouped_gen_tbl")

# Calculate Tajima's D across a window of 3 SNPs, with a step size of 2 SNPs
example_gt %>%
  windows_pop_tajimas_d(
    window_size = 3, step_size = 2,
    size_unit = "snp", min_loci = 2
  )
#> # A tibble: 3 × 7
#>   chromosome start   end n_loci  pop1   pop2    pop3
#>   <chr>      <dbl> <dbl>  <dbl> <dbl>  <dbl>   <dbl>
#> 1 chr1           1     3      3 1.03  -0.612  -0.710
#> 2 chr1           3     5      2 2.04  -0.612  -0.612
#> 3 chr2           1     3      2 0.311 -0.710 Inf    
```
