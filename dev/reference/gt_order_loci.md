# Order the loci table of a gen_tibble

This function reorders the loci table so that positions within a
chromosome are sequential. It also re-saves the genotypes into a new
file backed matrix with the new order, so that it can be used by
functions such as
[`loci_ld_clump()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/loci_ld_clump.md)
and
[`gt_pca_autoSVD()`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gt_pca_autoSVD.md).
If the loci table is already ordered, the original `gen_tibble` is
returned. This function will update the backingfiles of the `gen_tibble`
and return the `gen_tibble` object, use `<-` as per the example provided
to ensure that the names of the newly updated backingfiles are stored in
the `gen_tibble` object.

## Usage

``` r
gt_order_loci(
  .x,
  use_current_table = FALSE,
  ignore_genetic_dist = TRUE,
  quiet = FALSE,
  ...
)
```

## Arguments

- .x:

  a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)

- use_current_table:

  boolean, if FALSE (the default), the table will be reordered; if TRUE,
  then the current loci table, which might have been reordered manually,
  will be used, but only if the positions within each chromosome are
  sequential

- ignore_genetic_dist:

  boolean to ignore the genetic distance when checking. Note that, if
  `genetic_dist` are being ignored and they are not sorted, the function
  will set them to zero to avoid problems with other software.

- quiet:

  boolean to suppress information about the files

- ...:

  other arguments

## Value

A
[gen_tibble](https://evolecolgroup.github.io/tidypopgen/dev/reference/gen_tibble.md)

## Examples

``` r
example_gt <- load_example_gt("gen_tbl") %>% select_loci(c(1, 5, 2, 6, 4, 3))

# Loci are in the wrong order
show_loci(example_gt)
#> # A tibble: 6 × 7
#>   big_index name  chromosome position genetic_dist allele_ref allele_alt
#>       <int> <chr> <fct>         <int>        <dbl> <chr>      <chr>     
#> 1         1 rs1   chr1              3            0 A          T         
#> 2         5 rs5   chr2             23            0 C          G         
#> 3         2 rs2   chr1              5            0 T          C         
#> 4         6 rs6   chr2            456            0 T          A         
#> 5         4 rs4   chr1            343            0 G          C         
#> 6         3 rs3   chr1             65            0 C          NA        

# Reorder the loci, ignoring genetic distance
example_gt_ordered <- gt_order_loci(example_gt, ignore_genetic_dist = TRUE)
#> 
#> gen_backing files updated, now
#> using FBM RDS: /tmp/Rtmp4cRZgd/file25a34c82b055_v2.rds
#> with FBM backing file: /tmp/Rtmp4cRZgd/file25a34c82b055_v2.bk
#> make sure that you do NOT delete those files!
#> to reload the gen_tibble in another session, use:
#> gt_load('/tmp/Rtmp4cRZgd/file25a34c82b055_v2.gt')

# Loci are now in the correct order
show_loci(example_gt_ordered)
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
