# Create a Quality Control report for individuals

Return QC information to assess loci (Observed heterozygosity and
missingness).

## Usage

``` r
qc_report_indiv(.x, ...)

# S3 method for class 'tbl_df'
qc_report_indiv(.x, kings_threshold = NULL, ...)

# S3 method for class 'grouped_df'
qc_report_indiv(.x, kings_threshold = NULL, ...)
```

## Arguments

- .x:

  either a
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)
  object or a grouped
  [`gen_tibble`](https://evolecolgroup.github.io/tidypopgen/reference/gen_tibble.md)
  (as obtained by using
  [`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html))

- ...:

  further arguments to pass

- kings_threshold:

  an optional numeric giving a KING kinship coefficient, or one of:

  - "first": removing first degree relatives, equivalent to a kinship
    coefficient of 0.177 or more

  - "second": removing second degree relatives, equivalent to a kinship
    coefficient of 0.088 or more

## Value

If no kings_threshold is provided, a tibble with 2 elements: het_obs and
missingness. If kings_threshold is provided, a tibble with 4 elements:
het_obs, missingness, id and to_keep.

## Details

Providing the parameter kings_threshold will return two additional
columns, 'id' containing the ID of individuals, and 'to_keep' a logical
vector describing whether the individual should be removed to retain the
largest possible set of individuals with no relationships above the
threshold. The calculated pairwise KING relationship matrix is also
returned as an attribute of 'to_keep'. The kings_threshold parameter can
be either a numeric KING kinship coefficient or a string of either
"first" or "second", to remove any first degree or second degree
relationships from the dataset. This second option is similar to using
–unrelated –degree 1 or –unrelated –degree 2 in KING.

## Examples

``` r
# Create a gen_tibble of lobster genotypes
bed_file <-
  system.file("extdata", "lobster", "lobster.bed", package = "tidypopgen")
example_gt <- gen_tibble(bed_file,
  backingfile = tempfile("lobsters"),
  quiet = TRUE
)

# Get QC report for individuals
example_gt %>% qc_report_indiv()
#> # A tibble: 176 × 2
#>    het_obs missingness
#>      <dbl>       <dbl>
#>  1  0.0769      0.506 
#>  2  0.429       0.0253
#>  3  0           0.494 
#>  4  0.0930      0.456 
#>  5  0.0222      0.430 
#>  6  0.0833      0.544 
#>  7  0.0476      0.468 
#>  8  0.329       0.114 
#>  9  0.217       0.127 
#> 10  0.0857      0.557 
#> # ℹ 166 more rows

# Get QC report with kinship filtering
example_gt %>% qc_report_indiv(kings_threshold = "first")
#> # A tibble: 176 × 4
#>    het_obs missingness to_keep id   
#>      <dbl>       <dbl> <lgl>   <chr>
#>  1  0.0769      0.506  FALSE   Ale04
#>  2  0.429       0.0253 FALSE   Ale05
#>  3  0           0.494  FALSE   Ale06
#>  4  0.0930      0.456  FALSE   Ale08
#>  5  0.0222      0.430  FALSE   Ale13
#>  6  0.0833      0.544  FALSE   Ale15
#>  7  0.0476      0.468  FALSE   Ale16
#>  8  0.329       0.114  FALSE   Ale17
#>  9  0.217       0.127  FALSE   Ale18
#> 10  0.0857      0.557  FALSE   Ale19
#> # ℹ 166 more rows
```
