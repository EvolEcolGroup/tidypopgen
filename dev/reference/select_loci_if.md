# The `select_if` verb for `loci`

An equivalent to
[`dplyr::select_if()`](https://dplyr.tidyverse.org/reference/select_all.html)
that works on the `genotype` column of a `gen_tibble`. This function has
access to the genotypes (and thus can work on summary statistics to
select), but not the names of the loci (look at
[`select_loci()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/select_loci.md)
for that feature.

## Usage

``` r
select_loci_if(.data, .sel_logical)
```

## Arguments

- .data:

  a `gen_tibble`

- .sel_logical:

  a logical vector of length equal to the number of loci, or an
  expression that will tidy evaluate to such a vector. Only loci for
  which .sel_logical is TRUE will be selected; NA will be treated as
  FALSE.

## Value

a subset of the list of loci in the `gen_tibble`

## Details

Note that the `select_loci_if` verb does not modify the backing FBM
files, but rather it subsets the list of loci to be used stored in the
`gen_tibble`.

## See also

[`dplyr::select_if()`](https://dplyr.tidyverse.org/reference/select_all.html)

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

# Select loci by chromosome
example_gt_subset <- example_gt %>%
  select_loci_if(loci_chromosomes(genotypes) == "chr2")
show_loci(example_gt_subset)
#> # A tibble: 2 × 7
#>   big_index name  chromosome position genetic_dist allele_ref allele_alt
#>       <int> <chr> <fct>         <int>        <dbl> <chr>      <chr>     
#> 1         5 rs5   chr2             23            0 C          G         
#> 2         6 rs6   chr2            456            0 T          A         

# Select loci by a summary statistic
example_gt_subset <- example_gt %>%
  select_loci_if(loci_maf(genotypes) > 0.2)
show_loci(example_gt_subset)
#> # A tibble: 5 × 7
#>   big_index name  chromosome position genetic_dist allele_ref allele_alt
#>       <int> <chr> <fct>         <int>        <dbl> <chr>      <chr>     
#> 1         1 rs1   chr1              3            0 A          T         
#> 2         2 rs2   chr1              5            0 T          C         
#> 3         4 rs4   chr1            343            0 G          C         
#> 4         5 rs5   chr2             23            0 C          G         
#> 5         6 rs6   chr2            456            0 T          A         
```
