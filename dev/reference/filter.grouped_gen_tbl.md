# A filter method for grouped `gen_tibble` objects

A filter method for grouped `gen_tibble` objects

## Usage

``` r
# S3 method for class 'grouped_gen_tbl'
filter(..., deparse.level = 1)
```

## Arguments

- ...:

  a gen_tibble and a data.frame or tibble

- deparse.level:

  an integer controlling the construction of column names.

## Value

a grouped `gen_tibble`

## Examples

``` r
test_gt <- load_example_gt("grouped_gen_tbl")
test_gt %>% filter(id %in% c("a", "c"))
#> # A gen_tibble: 6 loci
#> # A tibble:     2 × 3
#> # Groups:       population [2]
#>   id    population  genotypes
#>   <chr> <chr>      <vctr_SNP>
#> 1 a     pop1        [1,1,...]
#> 2 c     pop2        [2,.,...]
test_gt <- load_example_gt("grouped_gen_tbl_sf")
test_gt %>% filter(id %in% c("a", "c"))
#> Simple feature collection with 2 features and 5 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 0 ymin: 49 xmax: 2 ymax: 51
#> Geodetic CRS:  WGS 84
#> # A gen_tibble: 6 loci
#> # A tibble:     2 × 6
#> # Groups:       population [2]
#>   id    population longitude latitude  genotypes    geometry
#> * <chr> <chr>          <dbl>    <dbl> <vctr_SNP> <POINT [°]>
#> 1 a     pop1               0       51  [1,1,...]      (0 51)
#> 2 c     pop2               2       49  [2,.,...]      (2 49)
```
