# Find duplicates in the loci table

This function finds duplicated SNPs by checking the positions within
each chromosome. It can return a list of duplicated SNPs or a logical
value indicating whether there are any duplicated loci.

## Usage

``` r
find_duplicated_loci(.x, error_on_false = FALSE, list_duplicates = TRUE, ...)
```

## Arguments

- .x:

  a vector of class `vctrs_bigSNP` (usually the `genotype` column of a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)
  object), or a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md).

- error_on_false:

  logical, if `TRUE` an error is thrown if duplicated loci are found.

- list_duplicates:

  logical, if `TRUE` returns duplicated SNP names.

- ...:

  other arguments passed to specific methods.

## Value

If `list_duplicates` is TRUE, returns a character vector of duplicated
loci names (character(0) when none). If `list_duplicates` is FALSE,
returns TRUE when no duplicates exist and FALSE when duplicates are
present. If `error_on_false` is TRUE and duplicates exist, an error is
thrown.

## Examples

``` r
example_gt <- load_example_gt("gen_tbl")
show_loci(example_gt) <- test_loci <- data.frame(
  big_index = c(1:6),
  name = paste0("rs", 1:6),
  chromosome = paste0("chr", c(1, 1, 1, 1, 1, 1)),
  position = as.integer(c(3, 3, 5, 65, 343, 46)),
  genetic_dist = as.double(rep(0, 6)),
  allele_ref = c("A", "T", "C", "G", "C", "T"),
  allele_alt = c("T", "C", NA, "C", "G", "A")
)

show_loci(example_gt)
#> # A tibble: 6 × 7
#>   big_index name  chromosome position genetic_dist allele_ref allele_alt
#>       <int> <chr> <fct>         <int>        <dbl> <chr>      <chr>     
#> 1         1 rs1   chr1              3            0 A          T         
#> 2         2 rs2   chr1              3            0 T          C         
#> 3         3 rs3   chr1              5            0 C          NA        
#> 4         4 rs4   chr1             65            0 G          C         
#> 5         5 rs5   chr1            343            0 C          G         
#> 6         6 rs6   chr1             46            0 T          A         

# Find which loci are duplicated
example_gt %>% find_duplicated_loci()
#> [1] "rs1" "rs2"
```
