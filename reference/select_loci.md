# The `select` verb for `loci`

An equivalent to
[`dplyr::select()`](https://dplyr.tidyverse.org/reference/select.html)
that works on the `genotype` column of a `gen_tibble`, using the
mini-grammar available for `tidyselect`. The `select`-like evaluation
only has access to the names of the loci (i.e. it can select only based
on names, not summary statistics of those loci; look at
[`select_loci_if()`](https://evolecolgroup.github.io/tidypopgen/reference/select_loci_if.md)
for that feature.

## Usage

``` r
select_loci(.data, .sel_arg)
```

## Arguments

- .data:

  a `gen_tibble`

- .sel_arg:

  one unquoted expression, using the mini-grammar of
  [`dplyr::select()`](https://dplyr.tidyverse.org/reference/select.html)
  to select loci. Variable names can be used as if they were positions
  in the data frame, so expressions like x:y can be used to select a
  range of variables.

## Value

a `gen_tibble` with a subset of the loci.

## Details

Note that the `select_loci` verb does not modify the backing FBM files,
but rather it subsets the list of loci to be used stored in the
`gen_tibble`.

## See also

[`dplyr::select()`](https://dplyr.tidyverse.org/reference/select.html)

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

# Select loci by name
example_gt_subset <- example_gt %>%
  select_loci(all_of(c("rs1", "rs2", "rs3")))
show_loci(example_gt_subset)
#> # A tibble: 3 × 7
#>   big_index name  chromosome position genetic_dist allele_ref allele_alt
#>       <int> <chr> <fct>         <int>        <dbl> <chr>      <chr>     
#> 1         1 rs1   chr1              3            0 A          T         
#> 2         2 rs2   chr1              5            0 T          C         
#> 3         3 rs3   chr1             65            0 C          NA        

# Select loci by index
example_gt_subset <- example_gt %>% select_loci(all_of(c(4, 2, 1)))
show_loci(example_gt_subset)
#> # A tibble: 3 × 7
#>   big_index name  chromosome position genetic_dist allele_ref allele_alt
#>       <int> <chr> <fct>         <int>        <dbl> <chr>      <chr>     
#> 1         4 rs4   chr1            343            0 G          C         
#> 2         2 rs2   chr1              5            0 T          C         
#> 3         1 rs1   chr1              3            0 A          T         
```
