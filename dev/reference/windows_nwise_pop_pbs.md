# Compute the Population Branch Statistics over a sliding window

The function computes the population branch statistics (PBS) for a
sliding window for each combination of populations at each locus. The
PBS is a measure of the genetic differentiation between one focal
population and two reference populations, and is used to identify
outlier loci that may be under selection.

## Usage

``` r
windows_nwise_pop_pbs(
  .x,
  type = c("matrix", "tidy"),
  fst_method = c("Hudson", "Nei87", "WC84"),
  return_fst = FALSE,
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
  followed by columns of pbs values for each population comparison.
  "tidy" returns a tidy tibble of the same data in 'long' format, where
  each row is one window for one population comparison.

- fst_method:

  the method to use for calculating Fst, one of 'Hudson', 'Nei87', and
  'WC84'. See
  [`pairwise_pop_fst()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pairwise_pop_fst.md)
  for details.

- return_fst:

  a logical value indicating whether to return the Fst values

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

- `pbs_a.b.c`: the PBS value for population a given b & c (there will be
  multiple such columns covering all 3 way combinations of populations
  in the grouped `gen_tibble` object)

- `fst_a.b`: the Fst value for population a and b, if `return_fst` is
  TRUE or a tidy tibble with the following columns:

- `chromosome`: the chromosome for the window

- `start`: the starting locus of the window

- `end`: the ending locus of the window

- `stat_name`: the name of populations used in the pbs calculation (e.g.
  "pbs_pop1.pop2.pop3"). If return_fst is TRUE, stat_name will also
  include "fst" calculations in the same column (e.g. "fst_pop1.pop2").

- `value`: the pbs value for the populations

## Examples

``` r
example_gt <- load_example_gt("grouped_gen_tbl")

# Calculate nwise pbs across a window of 3 SNPs, with a step size of 2 SNPs
example_gt %>%
  windows_nwise_pop_pbs(
    window_size = 3, step_size = 2,
    size_unit = "snp", min_loci = 2
  )
#> # A tibble: 3 × 9
#>   chromosome start   end pbs_pop1.pop2.pop3 pbs_pop2.pop1.pop3
#>   <chr>      <dbl> <dbl>              <dbl>              <dbl>
#> 1 chr1           1     3             0.0930             0.231 
#> 2 chr1           3     5             0.0486            -0.203 
#> 3 chr2           1     3            -0.0630            -0.0854
#> # ℹ 4 more variables: pbs_pop3.pop1.pop2 <dbl>, pbsn1_pop1.pop2.pop3 <dbl>,
#> #   pbsn1_pop2.pop1.pop3 <dbl>, pbsn1_pop3.pop1.pop2 <dbl>
```
