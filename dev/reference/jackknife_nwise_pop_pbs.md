# Compute the Population Branch Statistics with block jackknife standard errors

The function computes the genome-wide population branch statistics (PBS)
for each combination of populations, together with its standard error
estimated with a weighted block jackknife. The PBS is a measure of the
genetic differentiation between one focal population and two reference
populations.

## Usage

``` r
jackknife_nwise_pop_pbs(
  .x,
  type = c("matrix", "tidy"),
  fst_method = c("Hudson", "Nei87", "WC84"),
  return_fst = FALSE,
  block_size,
  size_unit = c("snp", "cM", "bp")
)
```

## Arguments

- .x:

  a grouped `gen_tibble` object

- type:

  type of object to return. One of "matrix" or "tidy". Default is
  "matrix". "matrix" returns a dataframe where each row is a quantity
  (estimate, jack_estimate, se, z and n_blocks), followed by columns for
  each population comparison. "tidy" returns a tidy tibble where each
  row is one population comparison.

- fst_method:

  the method to use for calculating Fst, one of 'Hudson', 'Nei87', and
  'WC84'. See
  [`pairwise_pop_fst()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pairwise_pop_fst.md)
  for details.

- return_fst:

  a logical value indicating whether to return the Fst values

- block_size:

  The size of the blocks to leave out.

- size_unit:

  Either "snp", "cM" or "bp". If "snp", the block size is in number of
  SNPs. If "cM", the block size is in centiMorgans. If "bp", the block
  size is in base pairs.

## Value

either a data frame with the following columns:

- `quantity`: the quantity reported in the row, one of `estimate`
  (computed on all loci), `jack_estimate` (the bias corrected jackknife
  estimate), `se` (the jackknife standard error), `z`
  (jack_estimate/se), and `n_blocks` (the number of blocks used)

- `fst_a.b`: the Fst value for population a and b, if `return_fst` is
  TRUE

- `pbs_a.b.c`: the PBS value for population a given b & c (there will be
  multiple such columns covering all 3 way combinations of populations
  in the grouped `gen_tibble` object)

- `pbsn1_a.b.c`: the normalised PBS value for population a given b & c
  or a tidy tibble with the following columns:

- `stat_name`: the name of populations used in the pbs calculation (e.g.
  "pbs_pop1.pop2.pop3"). If return_fst is TRUE, stat_name will also
  include "fst" calculations in the same column (e.g. "fst_pop1.pop2").

- `estimate`, `jack_estimate`, `se`, `z`, `n_blocks`: as described above

## Details

The genome is divided into non-overlapping blocks, which never span more
than one chromosome, and pairwise Fst (and thus PBS) is recomputed
leaving out one block at a time. The standard error is computed with the
weighted delete-m jackknife of Busing et al. (1999), where each block is
weighted by the number of loci it contains. See
[`jackknife_stats_generic()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/jackknife_stats_generic.md)
for details.

## References

Busing FMTA, Meijer E, Van Der Leeden R (1999) Delete-m jackknife for
unequal m. Statistics and Computing 9:3-8.

## Examples

``` r
example_gt <- load_example_gt("grouped_gen_tbl")

# Calculate nwise pbs with blocks of 2 SNPs
example_gt %>%
  jackknife_nwise_pop_pbs(
    block_size = 2,
    size_unit = "snp"
  )
#> # A tibble: 5 × 7
#>   quantity      pbs_pop1.pop2.pop3 pbs_pop2.pop1.pop3 pbs_pop3.pop1.pop2
#>   <chr>                      <dbl>              <dbl>              <dbl>
#> 1 estimate                  0.0461            -0.0110             0.0110
#> 2 jack_estimate             0.0743            -0.0274            -0.0512
#> 3 se                        0.0436             0.139              0.284 
#> 4 z                         1.70              -0.197             -0.180 
#> 5 n_blocks                  3                  3                  3     
#> # ℹ 3 more variables: pbsn1_pop1.pop2.pop3 <dbl>, pbsn1_pop2.pop1.pop3 <dbl>,
#> #   pbsn1_pop3.pop1.pop2 <dbl>
```
