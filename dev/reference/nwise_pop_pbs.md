# Compute the Population Branch Statistics for each combination of populations

The function computes the population branch statistics (PBS) for each
combination of populations at each locus. The PBS is a measure of the
genetic differentiation between one focal population and two reference
populations, and is used to identify outlier loci that may be under
selection.

## Usage

``` r
nwise_pop_pbs(
  .x,
  type = c("tidy", "matrix"),
  fst_method = c("Hudson", "Nei87", "WC84"),
  return_fst = FALSE
)
```

## Arguments

- .x:

  A grouped `gen_tibble`

- type:

  type of object to return. One of "tidy" or "matrix". Default is
  "tidy".

- fst_method:

  the method to use for calculating Fst, one of 'Hudson', 'Nei87', and
  'WC84'. See
  [`pairwise_pop_fst()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/pairwise_pop_fst.md)
  for details.

- return_fst:

  A logical value indicating whether to return the Fst values along with
  the PBS values. Default is `FALSE`.

## Value

Either a matrix with locus ID as rownames and the following columns:

- `pbs_a.b.c`: the PBS value for population a given b & c (there will be
  multiple such columns covering all 3 way combinations of populations
  in the grouped `gen_tibble` object)

- `pbsn1_a.b.c`: the normalized PBS value for population a given b & c.

- `fst_a.b`: the Fst value for population a and b, if `return_fst` is
  TRUE or a tidy tibble with the following columns:

- `loci`: the locus ID

- `stat_name`: the name of populations used in the pbs calculation (e.g.
  "pbs_pop1.pop2.pop3"). If return_fst is TRUE, stat_name will also
  include "fst" calculations in the same column (e.g. "fst_pop1.pop2").

- `value`: the pbs value for the populations

## References

Yi X, et al. (2010) Sequencing of 50 human exomes reveals adaptation to
high altitude. Science 329: 75-78.

## Examples

``` r
example_gt <- load_example_gt()

# We can compute the PBS for all populations using "Hudson" method
example_gt %>%
  group_by(population) %>%
  nwise_pop_pbs(fst_method = "Hudson")
#> # A tibble: 36 × 3
#>    loci  stat_name               value
#>    <chr> <chr>                   <dbl>
#>  1 rs1   pbs_pop1.pop2.pop3   -0.199  
#>  2 rs1   pbs_pop2.pop1.pop3   -0.0164 
#>  3 rs1   pbs_pop3.pop1.pop2    1.12   
#>  4 rs1   pbsn1_pop1.pop2.pop3 -0.105  
#>  5 rs1   pbsn1_pop2.pop1.pop3 -0.00863
#>  6 rs1   pbsn1_pop3.pop1.pop2  0.587  
#>  7 rs2   pbs_pop1.pop2.pop3    0.519  
#>  8 rs2   pbs_pop2.pop1.pop3    0.397  
#>  9 rs2   pbs_pop3.pop1.pop2   -0.397  
#> 10 rs2   pbsn1_pop1.pop2.pop3  0.342  
#> # ℹ 26 more rows
```
