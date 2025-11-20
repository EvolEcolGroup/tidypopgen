# Show the loci information of a `gen_tibble`

Extract and set the information on loci from a `gen_tibble`.

## Usage

``` r
show_loci(.x, ...)

# S3 method for class 'tbl_df'
show_loci(.x, ...)

# S3 method for class 'vctrs_bigSNP'
show_loci(.x, ...)

show_loci(.x) <- value

# S3 method for class 'tbl_df'
show_loci(.x) <- value

# S3 method for class 'vctrs_bigSNP'
show_loci(.x) <- value
```

## Arguments

- .x:

  a vector of class `vctrs_bigSNP` (usually the `genotype` column of a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)
  object), or a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md).

- ...:

  currently unused.

- value:

  a data.frame or tibble of loci information to replace the current one.

## Value

a [`tibble::tibble`](https://tibble.tidyverse.org/reference/tibble.html)
of information (see
[`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)
for details on compulsory columns that will always be present)

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")

example_gt %>% show_loci()
#> # A tibble: 6 × 7
#>   big_index name  chromosome position genetic_dist allele_ref allele_alt
#>       <int> <chr> <fct>         <int>        <dbl> <chr>      <chr>     
#> 1         1 rs1   chr1              3            0 A          T         
#> 2         2 rs2   chr1              5            0 T          C         
#> 3         3 rs3   chr1             65            0 C          NA        
#> 4         4 rs4   chr1            343            0 G          C         
#> 5         5 rs5   chr2             23            0 C          G         
#> 6         6 rs6   chr2            456            0 T          A         
```
