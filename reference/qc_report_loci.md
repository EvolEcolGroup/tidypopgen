# Create a Quality Control report for loci

Return QC information to assess loci (MAF, missingness and HWE test).
For pseudohaploid data, HWE test is not calculated.

## Usage

``` r
qc_report_loci(.x, ...)

# S3 method for class 'tbl_df'
qc_report_loci(.x, ...)

# S3 method for class 'grouped_df'
qc_report_loci(.x, ...)
```

## Arguments

- .x:

  a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)
  object.

- ...:

  currently unused

## Value

either a tibble with 3 elements (maf, missingness and hwe_p). For
pseudohaploid data, a tibble with 2 elements (maf and missingness).

## Examples

``` r
# Create a gen_tibble of lobster genotypes
bed_file <-
  system.file("extdata", "lobster", "lobster.bed", package = "tidypopgen")
example_gt <- gen_tibble(bed_file,
  backingfile = tempfile("lobsters"),
  quiet = TRUE
)

# Get a QC report for the loci
example_gt %>% qc_report_loci()
#> This gen_tibble is not grouped. For Hardy-Weinberg equilibrium, `qc_report_loci()` will assume individuals are part of the same population and HWE test p-values will be calculated across all individuals. If you wish to calculate HWE p-values within populations or groups, please use`group_by()` before calling `qc_report_loci()`.
#> # A tibble: 79 × 4
#>    snp_id    maf missingness       hwe_p
#>    <chr>   <dbl>       <dbl>       <dbl>
#>  1 rs3441  0.328      0.0227 0.440      
#>  2 rs4173  0.424      0.0341 0.584      
#>  3 rs6157  0.278      0.0284 0.00293    
#>  4 rs7502  0.423      0.0455 0.479      
#>  5 rs7892  0.101      0.0170 0.305      
#>  6 rs9441  0.307      0.0739 0.784      
#>  7 rs11071 0.123      0.0284 0.00299    
#>  8 rs11183 0.458      0.0568 0.102      
#>  9 rs11291 0.244      0.0455 0.000000944
#> 10 rs12971 0.312      0.0625 0.0242     
#> # ℹ 69 more rows

# Group by population to calculate HWE within populations
example_gt <- example_gt %>% group_by(population)
example_gt %>% qc_report_loci()
#> # A tibble: 79 × 4
#>    snp_id    maf missingness  hwe_p
#>    <chr>   <dbl>       <dbl>  <dbl>
#>  1 rs3441  0.328      0.0227 0.850 
#>  2 rs4173  0.424      0.0341 0.546 
#>  3 rs6157  0.278      0.0284 0.0834
#>  4 rs7502  0.423      0.0455 1.17  
#>  5 rs7892  0.101      0.0170 0.298 
#>  6 rs9441  0.307      0.0739 1.83  
#>  7 rs11071 0.123      0.0284 0.0143
#>  8 rs11183 0.458      0.0568 1.15  
#>  9 rs11291 0.244      0.0455 0.391 
#> 10 rs12971 0.312      0.0625 0.112 
#> # ℹ 69 more rows
```
