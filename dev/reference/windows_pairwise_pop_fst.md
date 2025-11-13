# Compute pairwise Fst for a sliding window

This function computes pairwise Fst for a sliding window across each
chromosome.

## Usage

``` r
windows_pairwise_pop_fst(
  .x,
  type = c("matrix", "tidy"),
  method = c("Hudson", "Nei87", "WC84"),
  window_size,
  step_size,
  size_unit = c("snp", "bp"),
  min_loci = 1,
  complete = FALSE
)
```

## Arguments

- .x:

  a grouped `gen_tibble` object

- type:

  type of object to return. One of "matrix" or "tidy". Default is
  "matrix". "matrix" returns a dataframe where each row is a window,
  followed by columns of Fst values for each pairwise population a and b
  comparison. "tidy" returns a tidy tibble of the same data in 'long'
  format, where each row is one window for one pairwise population a and
  b comparison.

- method:

  the method to use for calculating Fst, one of 'Hudson', 'Nei87', and
  'WC84'. See
  [`pairwise_pop_fst()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pairwise_pop_fst.md)
  for details.

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

either a data frame with the following columns:

- `chromosome`: the chromosome for the window

- `start`: the starting locus of the window

- `end`: the ending locus of the window

- `fst_a.b`: the pairwise Fst value for the population a and b (there
  will be multiple such columns if there are more than two populations)
  or a tidy tibble with the following columns:

- `chromosome`: the chromosome for the window

- `start`: the starting locus of the window

- `end`: the ending locus of the window

- `stat_name`: the name of population a and b used in the pairwise Fst
  calculation (e.g. "fst_pop1.pop2")

- `value`: the pairwise Fst value for the population a and b

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

example_gt %>%
  group_by(population) %>%
  windows_pairwise_pop_fst(
    window_size = 3, step_size = 2,
    size_unit = "snp", min_loci = 2
  )
#> # A tibble: 3 × 6
#>   chromosome start   end fst_pop1.pop2 fst_pop1.pop3 fst_pop2.pop3
#>   <chr>      <dbl> <dbl>         <dbl>         <dbl>         <dbl>
#> 1 chr1           1     3         0.277         0.311           0.4
#> 2 chr1           3     5        -0.167         0.222           0  
#> 3 chr2           1     3        -0.16         -0.467          -0.5
```
